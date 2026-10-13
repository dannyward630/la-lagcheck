import contextlib
import csv
from datetime import datetime, timedelta
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from examples.make_iolite_demo import make_demo
from lagcheck.cli import digest, main, run
from lagcheck import __version__
from lagcheck.iolite import export_iolite
from lagcheck.iolite_adapter import adapt
from lagcheck.models import InputError
from lagcheck.verification import verify_written_log



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

class IoliteWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.inputs = make_demo(self.root / 'inputs')
        self.canonical = self.root / 'canonical'
        code, self.report = run(self.inputs / 'config.json', self.canonical, True, 'affine',
                                '2026-01-01 00:00:00.000000', 'Known synthetic first-sample coordinate')
        self.assertEqual(code, 0, self.report.get('post_export_verification'))
        self.out = self.root / 'iolite'

    def test_cli_publishes_diagnosed_bytes_and_preserves_sources_without_claiming_gui(self):
        sources = list(self.inputs.glob('*')) + list(self.canonical.glob('*'))
        before = {str(path): digest(path) for path in sources if path.is_file()}
        console = io.StringIO()
        with contextlib.redirect_stdout(console):
            self.assertEqual(main(['iolite-export', str(self.canonical), '--output', str(self.out)]), 0)
        summary = json.loads(console.getvalue())
        self.assertEqual(summary['version'], __version__)
        self.assertFalse(summary['importer_roundtrip_tested'])
        self.assertFalse(summary['gui_workflow_tested'])
        final = Path(summary['usable_file'])
        post = json.loads((self.out / 'iolite_timing_verification.json').read_text(encoding='utf-8'))
        manifest = json.loads((self.out / 'correction_manifest.json').read_text(encoding='utf-8'))
        receipt = json.loads((self.out / 'format_stage/compatibility_receipt.json').read_text(encoding='utf-8'))
        self.assertEqual(final.read_bytes(), Path(receipt['output_path']).read_bytes())
        self.assertEqual(digest(final), post['diagnostic_input_sha256'])
        self.assertEqual(digest(final), manifest['output_sha256'])
        self.assertTrue(manifest['compatibility_adapter'])
        self.assertTrue(post['passed'])
        self.assertFalse(post['automatic_realignment_used'])
        self.assertEqual(post['alignment_offset_s'], 0.)
        self.assertFalse(receipt['rate_or_offset_applied'])
        self.assertTrue(receipt['non_time_cells_equal_by_source_record'])
        self.assertEqual(receipt['firing_interval_count'], 240)
        self.assertEqual(before, {path: digest(path) for path in before})

    def test_rejects_canonical_csv_changed_after_completed_manifest(self):
        path = self.canonical / 'corrected_laser_log.csv'
        path.write_bytes(path.read_bytes() + b'\n')
        with self.assertRaisesRegex(InputError, 'does not match'):
            export_iolite(self.canonical, self.out)
        self.assertFalse(self.out.exists())

    def test_rejects_incomplete_canonical_run(self):
        (self.canonical / 'INCOMPLETE.json').write_text('{}')
        with self.assertRaisesRegex(InputError, 'INCOMPLETE'):
            export_iolite(self.canonical, self.out)
        self.assertFalse(self.out.exists())

    def test_rejects_existing_destination_without_overwrite(self):
        self.out.mkdir()
        sentinel = self.out / 'keep.txt'
        sentinel.write_bytes(b'keep')
        with self.assertRaisesRegex(InputError, 'must not already exist'):
            export_iolite(self.canonical, self.out)
        self.assertEqual(sentinel.read_bytes(), b'keep')

    def test_rejects_changed_original_ms(self):
        signal = self.inputs / 'signal.csv'
        signal.write_bytes(signal.read_bytes() + b'\n')
        with self.assertRaisesRegex(InputError, 'ms input changed'):
            export_iolite(self.canonical, self.out)
        self.assertFalse(self.out.exists())

    def test_rejects_conflicting_canonical_time_mapping(self):
        path = self.canonical / 'correction_manifest.json'
        manifest = json.loads(path.read_text(encoding='utf-8'))
        manifest['export_plan']['coordinate_reference']['target_origin_seconds'] = '0'
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(InputError, 'fixed export mapping'):
            export_iolite(self.canonical, self.out)
        self.assertFalse(self.out.exists())

    def test_truthful_hash_does_not_hide_wrong_affine_constant(self):
        def shifted(*args, **kwargs):
            receipt = adapt(*args, **kwargs)
            candidate = Path(receipt['output_path'])
            with candidate.open(newline='', encoding='utf-8-sig') as stream:
                rows = list(csv.reader(stream))
            for row in rows[1:]:
                stamp = datetime.fromisoformat(row[0]) + timedelta(seconds=4)
                row[0] = stamp.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
            with candidate.open('w', newline='', encoding='utf-8') as stream:
                csv.writer(stream).writerows(rows)
            receipt['output_sha256'] = digest(candidate)
            return receipt
        with patch('lagcheck.iolite.adapt', side_effect=shifted):
            code, result = export_iolite(self.canonical, self.out)
        self.assertEqual(code, 2)
        post = result['post_export_verification']
        self.assertFalse(post['checks']['fixed_reference_alignment_within_tolerance'])
        self.assertEqual(post['alignment_offset_s'], 0.)
        self.assertFalse(post['automatic_realignment_used'])
        self.assertGreater(post['fixed_reference_max_lag_s'], 3.)
        self.assertIsNone(result['usable_file'])
        self.assertFalse((self.out / 'correction_manifest.json').exists())
        self.assertFalse((self.out / 'corrected_laser_log_iolite.csv').exists())

    def test_candidate_change_after_diagnostic_cannot_be_published(self):
        def changed(*args, **kwargs):
            result = verify_written_log(*args, **kwargs)
            path = Path(args[3])
            path.write_bytes(path.read_bytes() + b'\n')
            return result
        with patch('lagcheck.iolite.verify_written_log', side_effect=changed):
            with self.assertRaisesRegex(InputError, 'candidate changed'):
                export_iolite(self.canonical, self.out)
        self.assertTrue((self.out / 'INCOMPLETE.json').exists())
        self.assertFalse((self.out / 'correction_manifest.json').exists())
        self.assertFalse((self.out / 'corrected_laser_log_iolite.csv').exists())

    def test_input_change_during_diagnostic_cannot_be_published(self):
        def changed(*args, **kwargs):
            result = verify_written_log(*args, **kwargs)
            path = self.inputs / 'laser.csv'
            path.write_bytes(path.read_bytes() + b'\n')
            return result
        with patch('lagcheck.iolite.verify_written_log', side_effect=changed):
            with self.assertRaisesRegex(InputError, 'changed during'):
                export_iolite(self.canonical, self.out)
        self.assertTrue((self.out / 'INCOMPLETE.json').exists())
        self.assertFalse((self.out / 'correction_manifest.json').exists())

    def test_published_compatibility_hash_prevents_second_correction(self):
        code, result = export_iolite(self.canonical, self.out)
        self.assertEqual(code, 0)
        config = json.loads((self.canonical / 'config.json').read_text(encoding='utf-8'))
        config['laser']['path'] = result['usable_file']
        with self.assertRaisesRegex(InputError, 'earlier correction manifest'):
            run(config, self.root / 'again', True)
        self.assertFalse((self.root / 'again').exists())

    def test_rate_only_uses_existing_mapping_without_realigning(self):
        rate_canonical = self.root / 'rate_canonical'
        code, _ = run(self.inputs / 'config.json', rate_canonical, True)
        self.assertEqual(code, 0)
        code, result = export_iolite(rate_canonical, self.out)
        self.assertEqual(code, 0, result['post_export_verification'])
        self.assertEqual(result['canonical_mode'], 'rate_only')
        self.assertFalse(result['rate_or_offset_applied'])
        self.assertFalse(result['post_export_verification']['automatic_realignment_used'])
        self.assertAlmostEqual(result['post_export_verification']['alignment_offset_s'], 3., places=3)
        self.assertIn('continuous-trajectory transition times', result['receiving_software_instruction'])
        self.assertIn('firing-interval endpoints', result['receiving_software_instruction'])
        self.assertFalse(result['receiver_interval_trajectory_agreement_verified'])
        self.assertFalse(result['importer_roundtrip_tested'])
        self.assertEqual(result['acquisition_identity']['sequence_id'], '1')


if __name__ == '__main__':
    unittest.main()
