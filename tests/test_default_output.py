"""The public default must produce inspectable evidence without hiding failures."""
import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from lagcheck.cli import run, main, digest
from lagcheck.models import InputError
from lagcheck.laser import load_laser
from lagcheck.preview import render_preview_pair
from lagcheck.visualization import _grid, _map_lines
from lagcheck.verification import verify_written_log
from examples.make_iolite_demo import make_demo
import test_cli


class DefaultOutputTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=Path(__file__).parent)
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def case(self,name='inputs',step=0.,**settings):
        root=self.root/name; root.mkdir()
        cfg=test_cli.WorkflowTests.make_case(None,root,step)
        cfg['analysis'].setdefault('skip_record_edge_lines', False)
        cfg['analysis'].update(settings)
        (root/'config.json').write_text(json.dumps(cfg),encoding='utf-8')
        return root,cfg

    def read(self,path):
        return json.loads(path.read_text(encoding='utf-8'))

    def test_cli_defaults_cover_supported_failed_and_unestimated(self):
        for name,settings in [('supported',{}),('failed',{'step':2.}),('unestimated',{'peak_min':.999})]:
            root,cfg=self.case(name,**settings)
            out=self.root/(name+'_out')
            with contextlib.redirect_stdout(io.StringIO()) as console:
                code=main(['run',str(root/'config.json'),'--output',str(out)])
            self.assertIn(code,(0,2))
            self.assertEqual(json.loads(console.getvalue())['preview']['status'],'written')
            self.assertGreater((out/'preview_pair.png').stat().st_size,1000)
            meta=self.read(out/'preview_pair.json')
            self.assertEqual(meta['view'],'line_time')
            self.assertEqual(meta['pixel_pitch_unit'],'s')
            if name!='supported': self.assertTrue(meta['warnings'])
            if name=='unestimated': self.assertIsNone(meta['applied_rate_ppm'])

    def test_no_rate_uses_declared_alignment_and_reference_is_explicit(self):
        root,cfg=self.case(alignment_offset_s=7.,peak_min=.999)
        for ref in (None,170.9):
            out=self.root/('raw' if ref is None else 'reference')
            code,r=run(root/'config.json',out,reference_rate_ppm=ref)
            meta=self.read(out/'preview_pair.json')
            self.assertEqual(meta['alignment_offset_s'],7.)
            self.assertEqual(meta['applied_rate_ppm'],ref)
            self.assertEqual(meta['rate_source'],'user_reference_rate' if ref is not None else None)
            self.assertFalse(meta['rate_fitted'])

    def test_actual_published_mapping_overrides_baseline(self):
        root,cfg=self.case()
        out=self.root/'out'
        code,r=run(root/'config.json',out,True,'affine','0','Known synthetic coordinate')
        self.assertEqual(code,2)
        meta=self.read(out/'preview_pair.json')
        self.assertEqual(meta['rate_source'],'published_provisional_log')
        self.assertEqual(meta['mapping_provenance']['sha256'],digest(r['provisional_file']))
        actual=load_laser(r['provisional_file'],r['provisional_export_verification']['output_config'])
        shift=float(actual.table['origin'])
        self.assertTrue(np.allclose(meta['mapped_interval_starts_s'],[v['Start Time']+shift for v in actual.lines],rtol=0,atol=1e-9))
        self.assertTrue(np.allclose(meta['mapped_interval_ends_s'],[v['End Time']+shift for v in actual.lines],rtol=0,atol=1e-9))
        self.assertTrue(any('Failed checks' in w for w in meta['warnings']))

    def test_zero_duration_is_recorded_but_negative_duration_remains_error(self):
        root,cfg=self.case()
        path=root/'laser.csv'
        rows=list(csv.reader(path.open(newline='')))
        rows.insert(1,['0','0','zero'])
        with path.open('w',newline='') as f: csv.writer(f).writerows(rows)
        original=digest(path)
        code,r=run(root/'config.json',self.root/'zero',True)
        self.assertIn(code,(0,2))
        audit=self.read(self.root/'zero/input_processing.json')
        self.assertEqual(audit['excluded_count'],1)
        self.assertEqual(audit['excluded_intervals'][0]['start_record'],1)
        self.assertEqual(audit['retained_interval_count'],240)
        self.assertEqual(r['original_input_check']['status'],'requires_processing')
        self.assertFalse(r['correction_exported'])
        self.assertTrue((self.root/'zero/provisional_manifest.json').exists())
        self.assertEqual(digest(path),original)
        rows[1][1]='-1'
        with path.open('w',newline='') as f: csv.writer(f).writerows(rows)
        with self.assertRaises(InputError): run(root/'config.json',self.root/'negative')
        meta=self.read(self.root/'negative/preview_pair.json')
        self.assertEqual(meta['view'],'signal_time')
        self.assertTrue((self.root/'negative/INCOMPLETE.json').exists())

    def test_empty_off_fragment_is_not_a_second_sample_image(self):
        root=make_demo(self.root/'demo')
        path=root/'laser.csv'
        with path.open('a',encoding='utf-8',newline='') as f:
            w=csv.writer(f)
            w.writerow(['2026-01-01 05:00:00','Other','2','10','Off','meta'])
            w.writerow(['2026-01-01 06:00:00','Synthetic A','4','10','Off','meta'])
        code,r=run(root/'config.json',self.root/'segments',True)
        self.assertEqual(code,2)
        entries=self.read(self.root/'segments/segments.json')['segments']
        self.assertEqual(len(entries),2)
        self.assertEqual(entries[1]['status'],'no_firing_fragment')
        self.assertIsNone(entries[1]['preview'])
        self.assertTrue((self.root/'segments/segment_001/preview_pair.png').exists())
        self.assertFalse((self.root/'segments/segment_002/preview_pair.png').exists())
        self.assertFalse((self.root/'segments/correction_manifest.json').exists())
        self.assertEqual(entries[1]['origin_shift_s'],21600.)

    def test_no_overlap_gives_signal_view_and_large_grid_is_bounded(self):
        root,cfg=self.case(alignment_offset_s=1e9)
        rows=list(csv.reader((root/'laser.csv').open(newline='')))
        with (root/'laser.csv').open('w',newline='') as f: csv.writer(f).writerows(rows[:20])
        code,r=run(root/'config.json',self.root/'empty')
        meta=self.read(self.root/'empty/preview_pair.json')
        self.assertEqual(meta['view'],'signal_time')
        self.assertTrue(any('No MS samples' in v for v in meta['warnings']))
        geometry=dict(along_start=np.array([0.,0.]),along_end=np.array([1e9,1e9]),across=np.array([0.,1e9]),pitch=.01,row_pitch=.01)
        geom,grid=_grid(geometry)
        self.assertLessEqual(grid[1]*grid[3],1_000_000)
        self.assertGreater(geom['downsample_factor'],1)

    def test_candidate_byte_failure_keeps_exit_one_and_raw_view(self):
        root,cfg=self.case()
        def corrupt(*args,**kwargs):
            result=verify_written_log(*args,**kwargs)
            self.assertFalse(result[0]['passed'])
            p=Path(args[3]); p.write_bytes(p.read_bytes()+b'\n')
            return result
        out=self.root/'bad'
        with patch('lagcheck.cli.verify_written_log',side_effect=corrupt),contextlib.redirect_stderr(io.StringIO()):
            code=main(['run',str(root/'config.json'),'--output',str(out),'--export','--export-mode','affine','--ms-first-sample-time','0','--reference-note','Known synthetic coordinate'])
        self.assertEqual(code,1)
        self.assertTrue((out/'INCOMPLETE.json').exists())
        self.assertFalse((out/'provisional_manifest.json').exists())
        self.assertEqual(self.read(out/'preview_pair.json')['view'],'signal_time')

    def test_no_preview_opt_out_and_wizard_default(self):
        root,cfg=self.case()
        with contextlib.redirect_stdout(io.StringIO()):
            code=main(['run',str(root/'config.json'),'--output',str(self.root/'off'),'--no-preview'])
        self.assertFalse((self.root/'off/preview_pair.png').exists())
        with patch('lagcheck.wizard.collect',return_value=(root/'config.json',self.root/'wizard')),contextlib.redirect_stdout(io.StringIO()):
            code=main(['wizard'])
        self.assertTrue((self.root/'wizard/preview_pair.png').exists())

    def test_changed_published_bytes_are_not_drawn_as_correction(self):
        root,cfg=self.case()
        out=self.root/'out'
        code,r=run(root/'config.json',out,True,preview=False)
        p=Path(r['usable_file']); p.write_bytes(p.read_bytes()+b'\n')
        with self.assertRaisesRegex(InputError,'Published log bytes changed'):
            render_preview_pair(out)
        self.assertFalse((out/'preview_pair.png').exists())

    def test_unregistered_no_rate_never_calls_zero_a_known_origin(self):
        root,cfg=self.case()
        cfg['analysis'].pop('alignment_offset_s')
        (root/'config.json').write_text(json.dumps(cfg),encoding='utf-8')
        rows=list(csv.reader((root/'laser.csv').open(newline='')))
        with (root/'laser.csv').open('w',newline='') as f: csv.writer(f).writerows(rows[:20])
        code,r=run(root/'config.json',self.root/'raw')
        meta=self.read(self.root/'raw/preview_pair.json')
        self.assertIsNone(meta['alignment_offset_s'])
        self.assertIsNone(meta['applied_rate_ppm'])
        self.assertIn('unregistered',meta['alignment_source'])
        self.assertTrue(any('not a verified origin' in v for v in meta['warnings']))

    def test_invalid_fragment_does_not_silence_valid_fragment_or_error(self):
        root=make_demo(self.root/'demo')
        with (root/'laser.csv').open('a',encoding='utf-8',newline='') as f:
            w=csv.writer(f)
            w.writerow(['2026-01-01 05:00:00','Other','2','10','Off','meta'])
            w.writerow(['2026-01-01 06:00:00','Synthetic A','4','10','On','meta'])
        code,r=run(root/'config.json',self.root/'segments',True)
        self.assertEqual(code,1)
        self.assertEqual(r['status'],'segmented_with_errors')
        self.assertTrue((self.root/'segments/segment_001/preview_pair.png').exists())
        self.assertTrue((self.root/'segments/segment_002/INCOMPLETE.json').exists())
        self.assertEqual(self.read(self.root/'segments/segment_002/preview_pair.json')['view'],'signal_time')
        self.assertFalse((self.root/'segments/segment_002/correction_manifest.json').exists())

    def test_integrity_failure_during_default_preview_invalidates_run(self):
        root,cfg=self.case()
        from lagcheck.preview import render_preview_pair as real_render
        def altered(out,channel=None):
            report=self.read(Path(out)/'report.json')
            p=Path(report['usable_file']); p.write_bytes(p.read_bytes()+b'\n')
            return real_render(out,channel)
        with patch('lagcheck.preview.render_preview_pair',side_effect=altered),contextlib.redirect_stderr(io.StringIO()):
            code,r=run(root/'config.json',self.root/'out',True)
        self.assertEqual(code,1)
        self.assertTrue((self.root/'out/INCOMPLETE.json').exists())
        self.assertEqual(r['preview']['status'],'failed')

    def test_diagonal_stage_intervals_use_line_time_with_warning(self):
        from test_preview import PreviewTests
        case=PreviewTests(); case.root=self.root
        case.esl_case()
        path=self.root/'laser.csv'
        with path.open(newline='') as f: rows=list(csv.reader(f))
        rows[3][-1]='1.0'
        with path.open('w',newline='') as f: csv.writer(f).writerows(rows)
        code,r=run(self.root/'config.json',self.root/'out')
        meta=self.read(self.root/'out/preview_pair.json')
        self.assertEqual(meta['view'],'line_time')
        self.assertIn('diagonal or mixed-direction',meta['geometry_note'])

    def test_unnamed_sequence_ranges_do_not_repeat_full_input(self):
        from lagcheck.processing import inspect_segments
        path=self.root/'laser.csv'
        path.write_text('Timestamp,Laser State,Sequence Number\n0,Off,1\n1,On,\n2,Off,\n3,Off,2\n4,On,\n5,Off,\n',encoding='utf-8')
        cfg=dict(format='laser_state_csv',time_kind='numeric',time_unit='s')
        parts=inspect_segments(path,cfg)['segments']
        self.assertEqual(len(parts),2)
        for part in parts:
            config=dict(cfg,record_range=[part['start_record'],part['end_record']])
            self.assertEqual(len(inspect_segments(path,config)['segments']),1)
            loaded=load_laser(path,config)
            self.assertEqual(len(loaded.lines),1)
            self.assertEqual(loaded.metadata['acquisition_identity']['sequence_id'],part['sequence_id'])

    def test_missing_laser_keeps_error_but_readable_ms_can_be_inspected(self):
        root,cfg=self.case()
        cfg['laser']['path']='missing_laser.csv'
        (root/'config.json').write_text(json.dumps(cfg),encoding='utf-8')
        with contextlib.redirect_stderr(io.StringIO()):
            code=main(['run',str(root/'config.json'),'--output',str(self.root/'missing')])
        self.assertEqual(code,1)
        self.assertTrue((self.root/'missing/INCOMPLETE.json').exists())
        self.assertEqual(self.read(self.root/'missing/preview_pair.json')['view'],'signal_time')
        self.assertFalse((self.root/'missing/correction_manifest.json').exists())

    def test_trailing_nan_gets_partial_inspection_but_stays_strict_error(self):
        root,cfg=self.case()
        path=root/'signal.csv'
        with path.open(newline='') as f: rows=list(csv.reader(f))
        rows[-1][1]='NaN'
        with path.open('w',newline='') as f: csv.writer(f).writerows(rows)
        before=digest(path)
        with contextlib.redirect_stderr(io.StringIO()):
            code=main(['run',str(root/'config.json'),'--output',str(self.root/'partial'),'--export'])
        self.assertEqual(code,1)
        audit=self.read(self.root/'partial/input_processing.json')
        self.assertFalse(audit['analysis_cleaning_performed'])
        self.assertEqual(audit['invalid_record_count'],1)
        self.assertEqual(audit['invalid_signal_records'][0]['data_record'],len(rows)-1)
        self.assertTrue((self.root/'partial/INCOMPLETE.json').exists())
        self.assertTrue(self.read(self.root/'partial/preview_pair.json')['inspection_only'])
        self.assertFalse((self.root/'partial/correction_manifest.json').exists())
        self.assertEqual(digest(path),before)

    def test_new_inspection_paths_bind_the_bytes_actually_read(self):
        from lagcheck.models import IntegrityError
        from lagcheck.visualization import render_partial_signal
        from lagcheck.processing import inspect_segments as real_inspect
        root,cfg=self.case()
        for kind in ('ms','laser'):
            cfg[kind]['path']=str(root/cfg[kind]['path'])
        hashes={k:digest(cfg[k]['path']) for k in ('ms','laser')}
        read=Path.read_bytes
        target=Path(cfg['ms']['path']).resolve()
        def different_bytes(p):
            return read(p)+(b'\n' if p.resolve()==target else b'')
        with patch.object(Path,'read_bytes',new=different_bytes),self.assertRaises(IntegrityError):
            render_partial_signal(self.root/'partial',cfg,hashes,'strict input failed')
        def different_segmentation(*args,**kwargs):
            result=real_inspect(*args,**kwargs); result['source_sha256']='0'*64; return result
        with patch('lagcheck.processing.inspect_segments',side_effect=different_segmentation),self.assertRaises(IntegrityError):
            run(root/'config.json',self.root/'mismatched',preview=False)
        self.assertFalse((self.root/'mismatched').exists())
