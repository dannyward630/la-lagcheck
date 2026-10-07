"""Provisional (not verified) corrected logs, new in 0.2.0.

When the timing checks do not all pass, --export still writes a provisional
log with the estimated rate (or a user reference rate), marked as not verified.
The diagnosis and the verified export are unchanged.
"""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from examples.make_iolite_demo import make_demo
from lagcheck.cli import digest, main, run
from lagcheck.iolite import export_iolite
from lagcheck.laser import load_laser
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

class ProvisionalExportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=Path(__file__).parent)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def case(self, step=0., **analysis):
        config = test_cli.WorkflowTests.make_case(None, self.root, step)
        config['analysis'].update(analysis)
        (self.root / 'config.json').write_text(json.dumps(config), encoding='utf-8')
        return config

    def test_failed_checks_write_marked_provisional_log(self):
        config = self.case(step=2.)
        before = {n: digest(self.root / n) for n in ('signal.csv', 'laser.csv')}
        code, r = run(self.root / 'config.json', self.root / 'out', True)
        out = self.root / 'out'
        self.assertEqual(code, 2)
        self.assertNotEqual(r['status'], 'supported')
        self.assertFalse(r['correction_exported'])
        self.assertIsNone(r['usable_file'])
        self.assertFalse((out / 'corrected_laser_log.csv').exists())
        self.assertFalse((out / 'correction_manifest.json').exists())
        prov = out / 'provisional_laser_log.csv'
        manifest = json.loads((out / 'provisional_manifest.json').read_text(encoding='utf-8'))
        self.assertTrue(prov.is_file())
        self.assertTrue(manifest['provisional'])
        self.assertFalse(manifest['verified'])
        self.assertFalse(manifest['full_export_diagnostic_passed'])
        self.assertEqual(manifest['rate_source'], 'acquisition_estimate')
        self.assertEqual(manifest['output_sha256'], digest(prov))
        self.assertTrue(manifest['failed_checks'])
        self.assertAlmostEqual(manifest['rate_ppm'], r['rate_ppm'])
        self.assertIn('NOT verified', r['correction_note'])
        self.assertIn('Provisional (not verified) log', (out / 'report.md').read_text(encoding='utf-8'))
        # The written timestamps follow t' = t0 + (1 + r)(t - t0).
        export_config = {k: v for k, v in r['provisional_export_verification']['output_config'].items() if k != 'path'}
        original = load_laser(self.root / 'laser.csv', {k: v for k, v in config['laser'].items() if k != 'path'})
        written = load_laser(prov, export_config)
        t0 = original.lines[0]['Start Time']
        factor = 1 + manifest['rate_ppm'] * 1e-6
        for a, b in zip(original.lines, written.lines):
            for key in ('Start Time', 'End Time'):
                self.assertAlmostEqual(b[key] - written.lines[0]['Start Time'], factor * (a[key] - t0), delta=2e-6)
        for n, h in before.items():
            self.assertEqual(digest(self.root / n), h)

    def test_no_provisional_option_keeps_previous_behaviour(self):
        self.case(step=2.)
        console = io.StringIO()
        with contextlib.redirect_stdout(console):
            code = main(['run', str(self.root / 'config.json'), '--output', str(self.root / 'out'),
                         '--export', '--no-provisional'])
        self.assertEqual(code, 2)
        self.assertIsNone(json.loads(console.getvalue())['provisional_file'])
        self.assertFalse((self.root / 'out' / 'provisional_laser_log.csv').exists())
        self.assertFalse((self.root / 'out' / 'provisional_manifest.json').exists())

    def test_supported_run_writes_no_provisional_file(self):
        self.case()
        code, r = run(self.root / 'config.json', self.root / 'out', True)
        self.assertEqual(code, 0, r)
        self.assertTrue(r['correction_exported'])
        self.assertIsNone(r.get('provisional_file'))
        self.assertFalse((self.root / 'out' / 'provisional_manifest.json').exists())

    def test_reference_rate_used_only_when_no_rate_can_be_estimated(self):
        self.case(peak_min=0.999)
        code, r = run(self.root / 'config.json', self.root / 'out_none', True)
        self.assertEqual(code, 2)
        self.assertIsNone(r['rate_ppm'])
        self.assertIsNone(r.get('provisional_file'))
        self.assertIn('no rate could be estimated', r['provisional_note'])
        console = io.StringIO()
        with contextlib.redirect_stdout(console):
            code = main(['run', str(self.root / 'config.json'), '--output', str(self.root / 'out_ref'),
                         '--export', '--reference-rate-ppm', '170.9'])
        self.assertEqual(code, 2)
        manifest = json.loads((self.root / 'out_ref' / 'provisional_manifest.json').read_text(encoding='utf-8'))
        self.assertEqual(manifest['rate_source'], 'user_reference_rate')
        self.assertAlmostEqual(manifest['rate_ppm'], 170.9)
        self.assertEqual(manifest['mode'], 'rate_only')
        self.assertFalse(manifest['verified'])

    def test_reference_rate_must_be_numeric(self):
        config = self.case()
        config['export'] = {'reference_rate_ppm': 'fast'}
        (self.root / 'config.json').write_text(json.dumps(config), encoding='utf-8')
        with self.assertRaises(Exception):
            run(self.root / 'config.json', self.root / 'out', True)


class ProvisionalSafeguardTests(unittest.TestCase):
    """Regression tests for the R21 review: no second stretch, read-back stage, configuration errors, byte binding."""
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=Path(__file__).parent)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def case(self, step=0., export=None, **analysis):
        config = test_cli.WorkflowTests.make_case(None, self.root, step)
        config['analysis'].update(analysis)
        if export is not None:
            config['export'] = export
        (self.root / 'config.json').write_text(json.dumps(config), encoding='utf-8')
        return config

    def test_provisional_log_cannot_be_stretched_again(self):
        config = self.case(peak_min=0.999)
        code, r = run(self.root / 'config.json', self.root / 'first', True, reference_rate_ppm=170.9)
        self.assertEqual(code, 2)
        prov = Path(r['provisional_file'])
        config['laser']['path'] = str(prov)
        (self.root / 'again.json').write_text(json.dumps(config), encoding='utf-8')
        with contextlib.redirect_stderr(io.StringIO()):
            code = main(['run', str(self.root / 'again.json'), '--output', str(self.root / 'second'),
                         '--export', '--reference-rate-ppm', '170.9'])
        self.assertEqual(code, 1)
        self.assertFalse((self.root / 'second').exists())
        # Diagnosis without export remains available.
        code, _ = run(self.root / 'again.json', self.root / 'diagnosis_only', False)
        self.assertIn(code, (0, 2))

    def test_failed_readback_publishes_checked_candidate_as_provisional(self):
        import test_export_verification as tev
        case = tev.ExportVerificationTests.case
        case(None, self.root)
        code, r = run(self.root / 'config.json', self.root / 'out', True, 'affine', '0', 'Known synthetic coordinate')
        out = self.root / 'out'
        self.assertEqual(code, 2)
        self.assertEqual(r['input_diagnostic_status'], 'supported')
        self.assertEqual(r['status'], 'export_qc_failed')
        self.assertFalse(r['correction_exported'])
        self.assertIsNone(r['usable_file'])
        manifest = json.loads((out / 'provisional_manifest.json').read_text(encoding='utf-8'))
        self.assertEqual(manifest['failure_stage'], 'readback')
        self.assertFalse(manifest['verified'])
        self.assertEqual(manifest['output_sha256'], digest(out / 'provisional_laser_log.csv'))
        self.assertEqual(digest(out / 'provisional_laser_log.csv'), digest(out / 'candidate_laser_log.csv'))
        self.assertEqual(manifest['mode'], 'affine')
        code, r2 = run(self.root / 'config.json', self.root / 'off', True, 'affine', '0', 'Known synthetic coordinate', None, False)
        self.assertEqual(code, 2)
        self.assertFalse((self.root / 'off' / 'provisional_manifest.json').exists())

    def test_candidate_changed_after_readback_is_execution_error(self):
        from unittest.mock import patch
        import test_export_verification as tev
        from lagcheck.models import InputError
        from lagcheck.verification import verify_written_log

        def change_after_readback(ms, laser, report, candidate, receipt, plan, analysis):
            result = verify_written_log(ms, laser, report, candidate, receipt, plan, analysis)
            self.assertFalse(result[0]['passed'])    # the real read-back completed and did not pass
            candidate.write_bytes(candidate.read_bytes() + b'\n')
            return result
        tev.ExportVerificationTests.case(None, self.root)
        with patch('lagcheck.cli.verify_written_log', side_effect=change_after_readback):
            with self.assertRaises(InputError):
                run(self.root / 'config.json', self.root / 'out', True, 'affine', '0', 'Known synthetic coordinate')
        out = self.root / 'out'
        self.assertTrue((out / 'INCOMPLETE.json').exists())
        for name in ('provisional_laser_log.csv', 'provisional_manifest.json', 'corrected_laser_log.csv', 'correction_manifest.json'):
            self.assertFalse((out / name).exists(), name)
        # Through the command line the same fault is an execution error (exit code 1), not a failed timing check (2).
        with patch('lagcheck.cli.verify_written_log', side_effect=change_after_readback), contextlib.redirect_stderr(io.StringIO()):
            code = main(['run', str(self.root / 'config.json'), '--output', str(self.root / 'out_cli'), '--export',
                         '--export-mode', 'affine', '--ms-first-sample-time', '0', '--reference-note', 'Known synthetic coordinate'])
        self.assertEqual(code, 1)
        self.assertTrue((self.root / 'out_cli' / 'INCOMPLETE.json').exists())
        self.assertFalse((self.root / 'out_cli' / 'provisional_manifest.json').exists())

    def test_configuration_errors_are_input_errors(self):
        from lagcheck.models import InputError
        for export, name in (({'mode': 'affin'}, 'mode'), ({'typod_setting': 1}, 'key'),
                             ({'mode': 'affine', 'ms_first_sample_at_laser_time': 'not a time', 'reference_note': 'x'}, 'coordinate')):
            self.case(step=2., export=export)
            with self.assertRaises(InputError, msg=name):
                run(self.root / 'config.json', self.root / ('out_' + name), True)
            self.assertFalse((self.root / ('out_' + name)).exists())

    def test_reference_rate_with_affine_request_falls_back_with_note(self):
        self.case(peak_min=0.999, export={'mode': 'affine', 'ms_first_sample_at_laser_time': '0', 'reference_note': 'Synthetic'})
        code, r = run(self.root / 'config.json', self.root / 'out', True, reference_rate_ppm=170.9)
        self.assertEqual(code, 2)
        self.assertEqual(r['provisional_mode'], 'rate_only')
        self.assertIn('rate-only provisional log was written instead', r['provisional_note'])
        self.assertIn(r['provisional_readback']['status'], ('not_run', 'timing_checks_failed', 'passed'))

    def test_changed_bytes_after_writing_are_execution_errors(self):
        from unittest.mock import patch
        from lagcheck.laser import export_corrected
        from lagcheck.models import InputError

        def corrupt(laser, path, rate, **kwargs):
            receipt = export_corrected(laser, path, rate, **kwargs)
            path.write_bytes(path.read_bytes().replace(b'\n', b'\n ', 1))
            return receipt
        self.case(peak_min=0.999)
        with patch('lagcheck.cli.export_corrected', side_effect=corrupt):
            with self.assertRaises(InputError):
                run(self.root / 'config.json', self.root / 'out', True, reference_rate_ppm=170.9)
        self.assertFalse((self.root / 'out' / 'provisional_manifest.json').exists())
        self.assertTrue((self.root / 'out' / 'INCOMPLETE.json').exists())

    def test_existing_estimate_is_not_replaced_by_reference_rate(self):
        self.case(step=2.)
        code, r = run(self.root / 'config.json', self.root / 'out', True, reference_rate_ppm=999.)
        self.assertEqual(r['provisional_rate_source'], 'acquisition_estimate')
        self.assertAlmostEqual(r['provisional_rate_ppm'], r['rate_ppm'])


class ProvisionalIoliteExportTests(unittest.TestCase):
    def test_iolite_export_keeps_provisional_marking(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root = Path(folder)
            inputs = make_demo(root / 'inputs')
            config = json.loads((inputs / 'config.json').read_text(encoding='utf-8'))
            config['analysis']['min_coverage'] = 1.0      # full coverage is not reached, so one check fails; the rate stays estimable
            (inputs / 'config.json').write_text(json.dumps(config), encoding='utf-8')
            code, r = run(inputs / 'config.json', root / 'canonical', True, 'affine',
                          '2026-01-01 00:00:00.000000', 'Known synthetic first-sample coordinate')
            self.assertEqual(code, 2)
            self.assertFalse(r['correction_exported'])
            self.assertEqual(r['provisional_mode'], 'affine')
            code, result = export_iolite(root / 'canonical', root / 'iolite')
            self.assertEqual(code, 2)
            self.assertEqual(result['status'], 'provisional_not_verified')
            self.assertFalse(result['correction_exported'])
            self.assertIsNone(result['usable_file'])
            final = Path(result['provisional_file'])
            self.assertEqual(final.name, 'provisional_laser_log_iolite.csv')
            manifest = json.loads((root / 'iolite' / 'provisional_manifest.json').read_text(encoding='utf-8'))
            self.assertTrue(manifest['provisional'])
            self.assertFalse(manifest['verified'])
            self.assertEqual(manifest['output_sha256'], digest(final))
            self.assertFalse((root / 'iolite' / 'correction_manifest.json').exists())
            self.assertIn('PROVISIONAL, NOT VERIFIED', result['receiving_software_instruction'])
            # The iolite copy is refused as input to a further export too, so the rate cannot be applied twice.
            from lagcheck.models import InputError
            again = dict(config, laser=dict(config['laser'], path=str(final)),
                         ms=dict(config['ms'], path=str((inputs / config['ms']['path']).resolve())))
            (root / 'again.json').write_text(json.dumps(again), encoding='utf-8')
            with self.assertRaisesRegex(InputError, 'earlier correction manifest'):
                run(root / 'again.json', root / 'again', True)
            self.assertFalse((root / 'again').exists())


if __name__ == '__main__':
    unittest.main()
