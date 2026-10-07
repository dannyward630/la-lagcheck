"""Regression for repeated sample names across ESL acquisition sequences."""
import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from examples.make_iolite_demo import make_demo
from lagcheck.cli import main, run
from lagcheck.laser import load_laser, export_corrected
from lagcheck.models import InputError



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

class AcquisitionIdentityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=Path(__file__).parent)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def load(self, body, **options):
        p = self.root/'laser.csv'
        p.write_text('Timestamp,Laser State,Comment,Sequence Number\n'+body, encoding='utf-8')
        return load_laser(p, dict(time_kind='numeric', time_unit='s', **options))

    def test_same_sample_in_two_sequences_is_rejected(self):
        with self.assertRaisesRegex(InputError, 'Multiple Sequence Number'):
            self.load('0,Off,A,1\n1,On,,\n2,Off,,\n3,Off,A,2\n4,On,,\n5,Off,,\n', sample_selector='A')

    def test_single_sequence_inheritance_survives_export_with_receipt(self):
        data = self.load('0,Off,A,7\n1,On,,\n2,Off,,\n')
        self.assertEqual(data.metadata['acquisition_identity']['sequence_id'], '7')
        self.assertEqual(data.metadata['config']['expected_sequence'], '7')
        receipt = export_corrected(data, self.root/'corrected.csv', 100)
        self.assertEqual(receipt['acquisition_identity']['sequence_id'], '7')
        reloaded = load_laser(receipt['output_path'], receipt['output_config'])
        self.assertEqual(reloaded.metadata['acquisition_identity'], receipt['acquisition_identity'])

    def test_undefined_leading_sequence_is_rejected(self):
        with self.assertRaisesRegex(InputError, 'undefined Sequence Number'):
            self.load('0,Off,A,\n1,On,,1\n2,Off,,\n')

    def test_entirely_empty_sequence_is_rejected(self):
        with self.assertRaisesRegex(InputError, 'undefined Sequence Number'):
            self.load('0,Off,A,\n1,On,,\n2,Off,,\n')

    def test_expected_sequence_mismatch_is_rejected(self):
        with self.assertRaisesRegex(InputError, 'expected_sequence does not match'):
            self.load('0,Off,A,1\n1,On,,\n2,Off,,\n', expected_sequence='2')

    def test_other_sample_sequences_do_not_contaminate_selection(self):
        data = self.load('0,Off,A,1\n1,On,,\n2,Off,,\n3,Off,B,2\n4,On,,\n5,Off,,\n', sample_selector='B')
        self.assertEqual(data.metadata['acquisition_identity']['sequence_id'], '2')
        self.assertEqual(len(data.lines), 1)

    def test_generic_input_without_identity_has_explicit_limited_scope(self):
        p = self.root/'intervals.csv'
        p.write_text('start,end\n0,1\n2,3\n', encoding='utf-8')
        config = dict(format='interval_csv', time_kind='numeric', time_unit='s', start_column='start', end_column='end')
        data = load_laser(p, config)
        self.assertFalse(data.metadata['acquisition_identity']['sequence_verified'])
        with self.assertRaisesRegex(InputError, 'requires a Sequence Number'):
            load_laser(p, dict(config, expected_sequence='1'))

    def test_cli_splits_240_line_mixed_sequence_without_joint_export(self):
        inputs = make_demo(self.root/'inputs')
        p = inputs/'laser.csv'
        with p.open(encoding='utf-8', newline='') as f:
            rows = list(csv.reader(f))
        rows[242][1:3] = ['Synthetic A', '2']
        with p.open('w', encoding='utf-8', newline='') as f:
            csv.writer(f).writerows(rows)
        original = p.read_bytes()
        with contextlib.redirect_stdout(io.StringIO()):
            code = main(['run', str(inputs/'config.json'), '--output', str(self.root/'out'), '--export'])
        self.assertEqual(code, 2)
        index=json.loads((self.root/'out/segments.json').read_text(encoding='utf-8'))
        self.assertEqual([s['sequence_id'] for s in index['segments']], ['1','2'])
        self.assertFalse((self.root/'out/correction_manifest.json').exists())
        self.assertGreater(index['segments'][1]['declared_alignment_offset_s'],3.)
        self.assertEqual(p.read_bytes(), original)

    def test_canonical_config_manifest_and_warning_preserve_identity(self):
        inputs = make_demo(self.root/'inputs')
        code, report = run(inputs/'config.json', self.root/'out', True)
        self.assertEqual(code, 0, report)
        cfg = json.loads((self.root/'out/config.json').read_text(encoding='utf-8'))
        manifest = json.loads((self.root/'out/correction_manifest.json').read_text(encoding='utf-8'))
        self.assertEqual(cfg['laser']['expected_sequence'], '1')
        self.assertEqual(manifest['acquisition_identity']['sequence_id'], '1')
        self.assertIn('continuous-trajectory transition times', manifest['instruction'])
        self.assertFalse(manifest['receiver_interval_trajectory_agreement_verified'])
        self.assertFalse(manifest['importer_roundtrip_tested'])
