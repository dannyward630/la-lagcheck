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
from lagcheck.laser import export_corrected, load_laser
from lagcheck.models import InputError
from lagcheck.verification import verify_written_log
import test_cli



# Core regression tests opt out of the default rendering; visual behavior is tested separately.
from lagcheck.cli import run as _run_default_preview, main as _main_default_preview
def run(*args, **kwargs):
    if len(args) < 7:
        kwargs.setdefault('preview', False)
    elif args[6] is None:
        args = list(args); args[6] = False
    return _run_default_preview(*args, **kwargs)
def main(args=None):
    if args and args[0] == 'run' and '--preview' not in args and '--no-preview' not in args:
        args = list(args) + ['--no-preview']
    return _main_default_preview(args)

class ExportVerificationTests(unittest.TestCase):
    def case(self, root, rate=100.):
        config = test_cli.WorkflowTests().make_case(root)
        config['analysis']['skip_record_edge_lines'] = False
        (root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        if rate == 0.:
            # A no-rate control with grid-aligned irregular pulses isolates
            # recovery of the 3 s constant from the separate file origins.
            rng=np.random.default_rng(148);lines=[];cursor=5.
            for _ in range(240):
                length=round(float(rng.uniform(5,18)),1)
                lines.append((cursor,round(cursor+length,1)))
                cursor=round(cursor+length+float(rng.uniform(3,10)),1)
            time=np.arange(0,cursor+30,.1);signal=np.full(len(time),2.)
            for start,end in lines:
                signal[round((start+3)*10):round((end+3)*10)]=1000.
            with (root/'signal.csv').open('w', newline='') as f:
                writer=csv.writer(f);writer.writerow(['time_ms','element','not_signal'])
                writer.writerows((t*1000,y,999) for t,y in zip(time,signal))
            with (root/'laser.csv').open('w', newline='') as f:
                writer=csv.writer(f);writer.writerow(['start','end','x'])
                writer.writerows((start,end,i) for i,(start,end) in enumerate(lines))
        return config

    def test_affine_restores_known_numeric_coordinate_without_double_offset(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder); config=self.case(root,0.)
            code,r=run(root/'config.json',root/'out',True,'affine','0','Synthetic MS first sample is 0 s')
            self.assertEqual(code,0,r.get('post_export_verification'))
            self.assertTrue(r['correction_exported'])
            post=r['post_export_verification']
            self.assertEqual(post['alignment_offset_s'],0.)
            self.assertFalse(post['automatic_realignment_used'])
            self.assertLess(post['fixed_reference_max_lag_s'],.2)
            corrected=load_laser(r['usable_file'],r['export_verification']['output_config'])
            # Original first laser start is 5 s; true signal start is 8 s.
            self.assertAlmostEqual(float(corrected.table['intervals'][0][0]),8.,delta=.05)
            self.assertAlmostEqual(float(r['export_plan']['offset_seconds']),3.,delta=.05)
            config['laser']['path']=r['usable_file']
            config['ms']['path']=str(root/'signal.csv')
            with self.assertRaisesRegex(InputError,'earlier correction manifest'):
                run(config,root/'again',True)
            self.assertFalse((root/'again').exists())

    def test_wrong_rate_with_valid_csv_and_truthful_hash_is_rejected(self):
        def bad_export(laser,path,rate,**kwargs):
            return export_corrected(laser,path,0,**kwargs)
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder); self.case(root)
            with patch('lagcheck.cli.export_corrected',side_effect=bad_export):
                code,r=run(root/'config.json',root/'out',True)
            self.assertEqual(code,2)
            self.assertEqual(r['status'],'export_qc_failed')
            self.assertTrue(r['exported_log_diagnostic_rerun'])
            self.assertFalse(r['correction_exported'])
            self.assertIsNone(r['usable_file'])
            self.assertFalse(r['post_export_verification']['checks']['residual_rate_below_threshold'])
            # The rerun can itself support a NEW correction; that must not
            # conceal the fact that the actual exported file is still wrong.
            rr=json.loads((root/'out/post_export_report.json').read_text())
            self.assertEqual(rr['status'],'supported')
            self.assertFalse((root/'out/corrected_laser_log.csv').exists())
            # A written rate that differs from the plan is not released as a provisional log either.
            self.assertFalse((root/'out/provisional_laser_log.csv').exists())
            self.assertFalse((root/'out/provisional_manifest.json').exists())
            self.assertTrue(r['provisional_note'].startswith('No provisional log'))

    def test_affine_grid_quantization_can_honestly_prevent_export(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);self.case(root)
            code,r=run(root/'config.json',root/'out',True,'affine','0','Known synthetic coordinate')
            self.assertEqual(r['input_diagnostic_status'],'supported')
            self.assertEqual(code,2)
            self.assertFalse(r['correction_exported'])
            self.assertLess(r['post_export_verification']['fixed_reference_max_lag_s'],.2)
            self.assertFalse(r['post_export_verification']['checks']['residual_rate_below_threshold'])
            self.assertFalse((root/'out/corrected_laser_log.csv').exists())

    def test_wrong_affine_constant_is_not_hidden_by_realigning(self):
        def bad_export(laser,path,rate,**kwargs):
            kwargs['offset_seconds']=float(kwargs['offset_seconds'])+4.
            return export_corrected(laser,path,rate,**kwargs)
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);self.case(root)
            with patch('lagcheck.cli.export_corrected',side_effect=bad_export):
                code,r=run(root/'config.json',root/'out',True,'affine','0','Synthetic reference')
            self.assertEqual(code,2)
            self.assertFalse(r['correction_exported'])
            post=r['post_export_verification']
            self.assertFalse(post['checks']['fixed_reference_alignment_within_tolerance'])
            self.assertGreater(post['fixed_reference_max_lag_s'],3.)
            self.assertEqual(post['alignment_offset_s'],0.)
            self.assertFalse(post['automatic_realignment_used'])
            # A written constant that differs from the plan is not released as a provisional log either.
            self.assertFalse((root/'out/provisional_laser_log.csv').exists())
            self.assertFalse((root/'out/provisional_manifest.json').exists())
            self.assertTrue(r['provisional_note'].startswith('No provisional log'))

    def test_corrupt_on_disk_bytes_are_tool_error_not_qc_pass(self):
        def corrupt(laser,path,rate,**kwargs):
            receipt=export_corrected(laser,path,rate,**kwargs)
            path.write_bytes(path.read_bytes()+b'CORRUPTED\n')
            return receipt
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);self.case(root)
            with patch('lagcheck.cli.export_corrected',side_effect=corrupt),contextlib.redirect_stderr(io.StringIO()):
                code=main(['run',str(root/'config.json'),'--output',str(root/'out'),'--export'])
            self.assertEqual(code,1)
            self.assertTrue((root/'out/INCOMPLETE.json').exists())
            self.assertFalse((root/'out/corrected_laser_log.csv').exists())
            self.assertFalse((root/'out/correction_manifest.json').exists())

    def test_affine_does_not_guess_missing_reference(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);self.case(root)
            with self.assertRaisesRegex(InputError,'no clock origin is guessed'):
                run(root/'config.json',root/'out',True,'affine')
            self.assertFalse((root/'out/candidate_laser_log.csv').exists())

    def test_input_change_during_export_verification_never_publishes_manifest(self):
        def changed(*args,**kwargs):
            result=verify_written_log(*args,**kwargs)
            source=Path(args[1].metadata['source_path'])
            source.write_bytes(source.read_bytes()+b'\n')
            return result
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);self.case(root)
            with patch('lagcheck.cli.verify_written_log',side_effect=changed),contextlib.redirect_stderr(io.StringIO()):
                code=main(['run',str(root/'config.json'),'--output',str(root/'out'),'--export'])
            self.assertEqual(code,1)
            self.assertTrue((root/'out/INCOMPLETE.json').exists())
            self.assertFalse((root/'out/correction_manifest.json').exists())
            self.assertFalse((root/'out/report.json').exists())


if __name__=='__main__': unittest.main()
