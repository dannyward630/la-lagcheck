"""Strict diagnosis and independently gated provisional input-placement retries."""
import contextlib
import copy
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from lagcheck.cli import run, main, digest
from lagcheck.estimator import diagnose
from lagcheck.fallback import try_alignment_fallback, REQUIRED_CHECKS
from lagcheck.models import InputError
from examples.make_iolite_demo import make_demo
import test_cli


class AlignmentFallbackTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=Path(__file__).parent)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def case(self, declared=False):
        cfg = test_cli.WorkflowTests.make_case(None, self.root)
        if not declared:
            del cfg['analysis']['alignment_offset_s']
        (self.root/'config.json').write_text(json.dumps(cfg), encoding='utf-8')
        return cfg

    def json(self, path):
        return json.loads(path.read_text(encoding='utf-8'))

    def diagnostic(self, old_rate=False, reject=False, mutate=False):
        def implementation(ms, laser, analysis):
            report, blocks = diagnose(ms, laser, analysis)
            if analysis.get('alignment_offset_s') is None:
                report = copy.deepcopy(report)
                report.update(status='qc_failed' if old_rate else 'insufficient_evidence',
                              reason='qc_failed' if old_rate else 'initial_alignment_not_identifiable',
                              rate_ppm=37. if old_rate else None)
                report['alignment'] = dict(ok=bool(old_rate), source='automatic_candidate', offset_s=42. if old_rate else None)
                report['checks'] = {k:False for k in REQUIRED_CHECKS} if old_rate else {}
                self.strict = copy.deepcopy(report)
                return report, blocks
            self.candidate = copy.deepcopy(report)
            if reject:
                report['checks']['baseline_heldout'] = False
            if mutate:
                path = Path(laser.metadata['input_path']) if 'input_path' in laser.metadata else self.root/'laser.csv'
                path.write_bytes(path.read_bytes()+b'\n')
            return report, blocks
        return implementation

    def test_adopted_export_preserves_strict_and_uses_candidate_mapping(self):
        self.case()
        before = {n:digest(self.root/n) for n in ('signal.csv','laser.csv')}
        with patch('lagcheck.cli.diagnose', side_effect=self.diagnostic()):
            code,r = run(self.root/'config.json',self.root/'out',True,reference_rate_ppm=777.)
        self.assertEqual(code,2)
        for key in ('status','reason','rate_ppm','alignment','checks'):
            self.assertEqual(r[key],self.strict[key],key)
        self.assertEqual(r['input_diagnostic_status'],self.strict['status'])
        fallback=r['alignment_fallback']
        self.assertTrue(fallback['adopted'],fallback)
        self.assertEqual(r['provisional_rate_ppm'],self.candidate['rate_ppm'])
        self.assertEqual(r['provisional_export_plan']['initial_alignment_offset_s'],0.)
        manifest=self.json(self.root/'out/provisional_manifest.json')
        self.assertEqual(manifest['rate_source'],'input_position_fallback')
        self.assertFalse(manifest['verified'])
        self.assertEqual(manifest['diagnostic_status'],self.strict['status'])
        self.assertFalse(r['correction_exported'])
        self.assertFalse((self.root/'out/correction_manifest.json').exists())
        self.assertEqual(manifest['output_sha256'],digest(r['provisional_file']))
        meta=self.json(self.root/'out/preview_pair.json')
        self.assertTrue(meta['alignment_fallback']['adopted'])
        self.assertEqual(meta['mapping_provenance']['estimate_source'],'input_position_fallback')
        self.assertEqual(meta['alignment_offset_s'],0.)
        self.assertTrue(any('NOT VERIFIED' in w for w in meta['warnings']))
        for n,h in before.items(): self.assertEqual(digest(self.root/n),h)

    def test_real_cli_default_retry_and_opt_out_preserve_strict_fields(self):
        cfg=self.case()
        cfg['analysis']['coarse_window_s']=.01  # the real automatic window cannot identify a start
        (self.root/'config.json').write_text(json.dumps(cfg),encoding='utf-8')
        reports=[]
        for label,extra in [('default',[]),('off',['--no-alignment-fallback'])]:
            command=[sys.executable,'-B','-X','utf8','-m','lagcheck','run',str(self.root/'config.json'),
                     '--output',str(self.root/label),'--export']+extra
            completed=subprocess.run(command,cwd=Path(__file__).resolve().parents[1],capture_output=True,text=True)
            self.assertEqual(completed.returncode,2,completed.stderr)
            reports.append(self.json(self.root/label/'report.json'))
            self.assertTrue((self.root/label/'preview_pair.png').exists())
        for key in ('status','reason','rate_ppm','alignment','checks','input_diagnostic_status'):
            self.assertEqual(reports[0].get(key),reports[1].get(key),key)
        self.assertEqual(reports[0]['reason'],'initial_alignment_not_identifiable')
        self.assertTrue(reports[0]['alignment_fallback']['adopted'])
        self.assertFalse(reports[1]['alignment_fallback']['attempted'])
        self.assertTrue((self.root/'default/provisional_manifest.json').exists())
        self.assertFalse((self.root/'off/provisional_manifest.json').exists())

    def test_bad_automatic_rate_is_replaced_only_in_provisional_outputs(self):
        self.case()
        with patch('lagcheck.cli.diagnose',side_effect=self.diagnostic(old_rate=True)):
            code,r=run(self.root/'config.json',self.root/'out',True)
        self.assertEqual(code,2)
        self.assertEqual(r['rate_ppm'],37.)
        self.assertEqual(r['alignment']['offset_s'],42.)
        self.assertEqual(r['alignment_fallback']['triggered_reason'],'automatic_estimate_failed_agreement_or_prediction')
        self.assertEqual(r['provisional_rate_ppm'],self.candidate['rate_ppm'])
        self.assertEqual(r['provisional_export_plan']['initial_alignment_offset_s'],0.)
        self.assertEqual(self.json(self.root/'out/preview_pair.json')['alignment_offset_s'],0.)

    def test_rejected_fallback_keeps_raw_view_and_records_why(self):
        self.case()
        with patch('lagcheck.cli.diagnose',side_effect=self.diagnostic(reject=True)):
            code,r=run(self.root/'config.json',self.root/'out',True)
        self.assertEqual(code,2)
        self.assertFalse(r['alignment_fallback']['adopted'])
        self.assertIn('baseline_heldout',r['alignment_fallback']['failed_adoption_checks'])
        self.assertFalse((self.root/'out/provisional_manifest.json').exists())
        meta=self.json(self.root/'out/preview_pair.json')
        self.assertIsNone(meta['applied_rate_ppm'])
        self.assertTrue(any('attempted but not adopted' in w for w in meta['warnings']))

    def test_no_export_still_estimates_and_draws_but_no_log(self):
        self.case()
        with patch('lagcheck.cli.diagnose',side_effect=self.diagnostic()):
            code,r=run(self.root/'config.json',self.root/'out')
        self.assertEqual(code,2)
        self.assertTrue(r['alignment_fallback']['adopted'])
        self.assertFalse(list((self.root/'out').glob('*laser_log.csv')))
        meta=self.json(self.root/'out/preview_pair.json')
        self.assertEqual(meta['rate_source'],'input_position_fallback')
        self.assertEqual(meta['applied_rate_ppm'],self.candidate['rate_ppm'])

    def test_explicit_alignment_and_disabling_switches_prevent_attempt(self):
        self.case(declared=True)
        code,r=run(self.root/'config.json',self.root/'declared',preview=False)
        self.assertFalse(r['alignment_fallback']['attempted'])
        cfg=self.json(self.root/'config.json'); del cfg['analysis']['alignment_offset_s']
        (self.root/'config.json').write_text(json.dumps(cfg),encoding='utf-8')
        for index,flag in enumerate(('--no-alignment-fallback','--no-provisional')):
            with patch('lagcheck.cli.diagnose',side_effect=self.diagnostic()) as call, contextlib.redirect_stdout(io.StringIO()):
                code=main(['run',str(self.root/'config.json'),'--output',str(self.root/str(index)),'--export','--no-preview',flag])
            self.assertEqual(code,2)
            self.assertEqual(call.call_count,1)
            self.assertFalse(self.json(self.root/str(index)/'report.json')['alignment_fallback']['attempted'])

    def test_no_preview_keeps_adopted_fallback_and_readback_never_upgrades(self):
        self.case()
        with patch('lagcheck.cli.diagnose',side_effect=self.diagnostic()):
            code,r=run(self.root/'config.json',self.root/'out',True,preview=False)
        self.assertEqual(code,2)
        self.assertTrue(r['alignment_fallback']['adopted'])
        self.assertEqual(r['provisional_readback']['status'],'passed')
        self.assertFalse((self.root/'out/preview_pair.png').exists())
        self.assertFalse(self.json(self.root/'out/provisional_manifest.json')['verified'])

    def test_rejected_retry_retains_existing_strict_rate_and_mapping(self):
        self.case()
        with patch('lagcheck.cli.diagnose',side_effect=self.diagnostic(old_rate=True,reject=True)):
            code,r=run(self.root/'config.json',self.root/'out',True)
        self.assertEqual(code,2)
        self.assertFalse(r['alignment_fallback']['adopted'])
        self.assertEqual(r['provisional_rate_ppm'],37.)
        self.assertEqual(r['provisional_rate_source'],'acquisition_estimate')
        self.assertEqual(r['provisional_export_plan']['initial_alignment_offset_s'],42.)
        meta=self.json(self.root/'out/preview_pair.json')
        self.assertEqual(meta['alignment_offset_s'],42.)
        self.assertEqual(meta['applied_rate_ppm'],37.)
        self.assertTrue(any('attempted but not adopted' in w for w in meta['warnings']))

    def test_iolite_export_uses_adopted_candidate_and_keeps_provisional_status(self):
        from lagcheck.iolite import export_iolite
        from lagcheck.verification import verify_written_log
        root=make_demo(self.root/'demo')
        cfg=self.json(root/'config.json'); del cfg['analysis']['alignment_offset_s']
        (root/'config.json').write_text(json.dumps(cfg),encoding='utf-8')
        for old_rate in (False,True):
            canonical=self.root/('canonical_'+str(old_rate))
            adapted=self.root/('adapted_'+str(old_rate))
            with patch('lagcheck.cli.diagnose',side_effect=self.diagnostic(old_rate=old_rate)):
                code,r=run(root/'config.json',canonical,True,preview=False)
            self.assertEqual(code,2)
            self.assertTrue(r['alignment_fallback']['adopted'])
            with patch('lagcheck.iolite.verify_written_log',wraps=verify_written_log) as actual_readback:
                code,result=export_iolite(canonical,adapted)
            self.assertEqual(code,2)
            actual_readback.assert_called_once()
            timing_report=actual_readback.call_args.args[2]
            self.assertEqual(timing_report,r['alignment_fallback']['diagnostic'])
            self.assertEqual(timing_report['alignment']['offset_s'],0.)
            self.assertNotEqual(r['alignment']['offset_s'],0.)
            self.assertFalse(result['correction_exported'])
            self.assertIsNotNone(result['post_export_verification'])
            manifest=self.json(adapted/'provisional_manifest.json')
            self.assertFalse(manifest['verified'])
            self.assertEqual(manifest['rate_source'],'input_position_fallback')
            self.assertEqual(manifest['alignment_fallback'],r['alignment_fallback'])
            self.assertEqual(manifest['output_sha256'],digest(result['provisional_file']))
            self.assertEqual(result['post_export_verification']['diagnostic_input_sha256'],manifest['output_sha256'])

    def test_each_required_gate_missing_or_false_and_nonfinite_rate_reject(self):
        strict=dict(reason='initial_alignment_not_identifiable',rate_ppm=None)
        good=dict(rate_ppm=170.,checks={k:True for k in REQUIRED_CHECKS},closure={'converged':False})
        accepted,_=try_alignment_fallback(None,None,strict,{},lambda *a:(good,[]))
        self.assertTrue(accepted['adopted'])  # convergence is not an extra gate
        cases=[]
        for key in REQUIRED_CHECKS:
            for value in (False,None,'missing'):
                bad=copy.deepcopy(good)
                if value=='missing': del bad['checks'][key]
                else: bad['checks'][key]=value
                cases.append((key+' '+str(value),bad))
        for rate in (None,float('nan'),float('inf'),-1000000.,True):
            cases.append(('rate '+str(rate),dict(good,rate_ppm=rate)))
        for label,candidate in cases:
            with self.subTest(case=label):
                result,_=try_alignment_fallback(None,None,strict,{},lambda *a:(candidate,[]))
                self.assertFalse(result['adopted'])

    def test_other_failure_does_not_trigger_and_retry_errors_propagate(self):
        from unittest.mock import Mock
        for reason,rate,checks in [('time_grid_invalid',None,{}),('coverage',170.,{'coverage':False}),
                                   ('sensitivity',170.,{'sensitivity_envelope':False})]:
            fn=Mock()
            r,_=try_alignment_fallback(None,None,dict(reason=reason,rate_ppm=rate,checks=checks),{},fn)
            self.assertFalse(r['attempted']); fn.assert_not_called()
        with self.assertRaises(InputError):
            try_alignment_fallback(None,None,dict(reason='initial_alignment_not_identifiable',rate_ppm=None),{},
                                   Mock(side_effect=InputError('invalid candidate input')))

    def test_input_mutation_during_retry_prevents_export(self):
        self.case()
        with patch('lagcheck.cli.diagnose',side_effect=self.diagnostic(mutate=True)):
            with self.assertRaisesRegex(InputError,'Input changed'):
                run(self.root/'config.json',self.root/'out',True,preview=False)
        self.assertFalse((self.root/'out/provisional_manifest.json').exists())

    def test_zero_duration_and_segment_children_use_retry_without_resetting_origin(self):
        root=make_demo(self.root/'demo')
        cfg=self.json(root/'config.json'); del cfg['analysis']['alignment_offset_s']
        (root/'config.json').write_text(json.dumps(cfg),encoding='utf-8')
        path=root/'laser.csv'
        rows=list(csv.reader(path.open(encoding='utf-8',newline='')))
        rows[1:1]=[['2026-01-01 00:00:00','Synthetic A','1','10','On','zero'],
                   ['2026-01-01 00:00:00','','','','Off','zero']]
        # A second selected fragment begins later on the same original MS axis.
        rows += [['2026-01-01 02:00:00','Other','2','10','Off','metadata'],
                 ['2026-01-01 03:00:00','Synthetic A','3','10','On','later'],
                 ['2026-01-01 03:00:01','','','','Off','later']]
        with path.open('w',encoding='utf-8',newline='') as f: csv.writer(f).writerows(rows)
        with patch('lagcheck.cli.diagnose',side_effect=self.diagnostic()):
            code,r=run(root/'config.json',self.root/'segments',True,preview=False)
        self.assertEqual(code,2)
        first=self.json(self.root/'segments/segment_001/report.json')
        later=self.json(self.root/'segments/segment_002/report.json')
        self.assertTrue(first['alignment_fallback']['adopted'])
        self.assertIn('zero_duration_intervals',first['original_input_check']['reasons'])
        self.assertTrue(later['alignment_fallback']['attempted'])
        self.assertEqual(later['alignment_fallback']['alignment_offset_s'],10800.)
        self.assertFalse(later['alignment_fallback']['adopted'])
        self.assertFalse(first['correction_exported'])


if __name__=='__main__': unittest.main()
