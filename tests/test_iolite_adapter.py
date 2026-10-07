"""Behavioral checks, including the two exact R005 canonical outputs."""
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys
import os
import tempfile
import unittest
from unittest.mock import patch

from lagcheck import iolite_adapter as adapter

HERE = Path(__file__).resolve().parent
SOFTWARE = HERE.parent



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

def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.reader(stream))


class AdapterTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="adapter_test_", dir=HERE)
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def case(self, rows=None):
        source = self.root / "canonical.csv"
        if rows is None:
            rows = [
                ["2026-01-01 00:00:00.000000", "A", "1", "10", "Off", " 01 "],
                ["2026-01-01 00:00:04.000000", "", "", "10", "Off", "planned"],
                ["2026-01-01 00:00:01.000500", "", "", "", "On", "first edge"],
                ["2026-01-01 00:00:02.001500", "", "", "", "Off", "last edge"],
            ]
        with source.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["Timestamp", " Comment", " Sequence Number", " Scan Velocity (um/s)", " Laser State", "X"])
            writer.writerows(rows)
        manifest = self.root / "source_manifest.json"
        manifest.write_text(json.dumps(dict(version="0.1.0", mode="affine",
                            output_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                            full_export_diagnostic_passed=True,
                            export_plan={"offset_seconds": "3.0"})), encoding="utf-8")
        return source, manifest

    def test_known_dates_stable_reordering_and_no_second_correction(self):
        source, manifest = self.case()
        before = source.read_bytes()
        result = adapter.adapt(source, manifest, self.root / "out")
        output = Path(result["output_path"])
        rows = read_csv(output)[1:]
        self.assertEqual([r[0] for r in rows], ["2026-01-01 00:00:00.000", "2026-01-01 00:00:01.000",
                                               "2026-01-01 00:00:02.002", "2026-01-01 00:00:04.000"])
        self.assertEqual([r[-1] for r in rows], [" 01 ", "first edge", "last edge", "planned"])
        self.assertEqual(result["source_record_by_output_record"], [1, 3, 4, 2])
        self.assertEqual(result["firing_endpoint_source_records"], [[3, 4]])
        self.assertFalse(result["rate_or_offset_applied"])
        self.assertFalse(result["full_export_diagnostic_passed"])
        self.assertFalse(result["importer_roundtrip_tested"])
        self.assertIsNone(result["usable_file"])
        self.assertEqual(source.read_bytes(), before)
        sys.path.insert(0, str(SOFTWARE))
        from lagcheck.cli import _reject_known_corrected_input
        from lagcheck.models import InputError
        with self.assertRaisesRegex(InputError, "earlier correction manifest"):
            _reject_known_corrected_input(output)
        with self.assertRaisesRegex(adapter.CompatibilityError, "already adapted"):
            adapter.adapt(output, output.parent / "correction_manifest.json", self.root / "again")

    def test_millisecond_rounding_crosses_midnight(self):
        source, manifest = self.case([
            ["2026-01-01 23:59:58.000000", "A", "1", "10", "Off", ""],
            ["2026-01-01 23:59:58.999500", "", "", "", "On", ""],
            ["2026-01-01 23:59:59.999500", "", "", "", "Off", ""],
        ])
        result = adapter.adapt(source, manifest, self.root / "out")
        rows = read_csv(result["output_path"])[1:]
        self.assertEqual(rows[-1][0], "2026-01-02 00:00:00.000")
        self.assertEqual(rows[-2][0], "2026-01-01 23:59:59.000")

    def test_multiple_samples_and_sequences_are_rejected(self):
        for col in (1, 2):
            with self.subTest(column=col):
                folder = self.root / str(col)
                folder.mkdir()
                source, manifest = self.case_in(folder, change=(col, "B" if col == 1 else "2"))
                with self.assertRaisesRegex(adapter.CompatibilityError, "Only one selected"):
                    adapter.adapt(source, manifest, folder / "out")
                self.assertFalse((folder / "out").exists())

    def case_in(self, folder, change=None):
        old = self.root
        self.root = folder
        try:
            source, manifest = self.case()
            rows = read_csv(source)
            if change:
                rows[-1][change[0]] = change[1]
                self.rewrite_source(source, manifest, rows)
            return source, manifest
        finally:
            self.root = old

    def rewrite_source(self, source, manifest, rows):
        with source.open("w", encoding="utf-8", newline="") as stream:
            csv.writer(stream).writerows(rows)
        data = json.loads(manifest.read_text())
        data["output_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
        manifest.write_text(json.dumps(data), encoding="utf-8")

    def test_sorting_that_closes_on_with_different_off_is_rejected(self):
        source, manifest = self.case()
        rows = read_csv(source)
        rows[2][0] = "2026-01-01 00:00:01.500000"
        self.rewrite_source(source, manifest, rows)
        with self.assertRaisesRegex(adapter.CompatibilityError, "paired firing endpoints"):
            adapter.adapt(source, manifest, self.root / "out")
        self.assertFalse((self.root / "out").exists())

    def test_firing_velocity_change_is_rejected(self):
        source, manifest = self.case()
        rows = read_csv(source)
        rows[2][3] = "20"
        self.rewrite_source(source, manifest, rows)
        with self.assertRaisesRegex(adapter.CompatibilityError, "scan velocity"):
            adapter.adapt(source, manifest, self.root / "out")

    def test_empty_sample_inheritance_is_materialized_before_sorting(self):
        source, manifest = self.case()
        rows = read_csv(source)
        rows[2][0] = "2025-12-31 23:59:59.000000"
        self.rewrite_source(source, manifest, rows)
        result = adapter.adapt(source, manifest, self.root / "out")
        self.assertTrue(result['sample_and_sequence_inheritance_equal_all_records'])
        self.assertFalse(result['non_time_cells_equal_by_source_record'])
        self.assertTrue(result['non_membership_cells_equal_by_source_record'])
        self.assertEqual(read_csv(result['output_path'])[1][1:3], ['A', '1'])
        self.assertEqual(result['source_record_by_output_record'], [2, 1, 3, 4])

    def test_precision_collapse_is_rejected(self):
        source, manifest = self.case([
            ["2026-01-01 00:00:00.000000", "A", "1", "10", "Off", ""],
            ["2026-01-01 00:00:01.000100", "", "", "", "On", ""],
            ["2026-01-01 00:00:01.000400", "", "", "", "Off", ""],
        ])
        with self.assertRaisesRegex(adapter.CompatibilityError, "nonpositive"):
            adapter.adapt(source, manifest, self.root / "out")
        self.assertFalse((self.root / "out").exists())

    def test_input_hash_mismatch_is_rejected(self):
        source, manifest = self.case()
        source.write_bytes(source.read_bytes() + b"\n")
        with self.assertRaisesRegex(adapter.CompatibilityError, "manifest hash"):
            adapter.adapt(source, manifest, self.root / "out")
        self.assertFalse((self.root / "out").exists())

    def test_source_changed_during_preparation_is_rejected(self):
        source, manifest = self.case()
        original = adapter.prepare
        def changed(payload):
            result = original(payload)
            source.write_bytes(source.read_bytes() + b"\n")
            return result
        with patch.object(adapter, "prepare", side_effect=changed):
            with self.assertRaisesRegex(adapter.CompatibilityError, "changed during preparation"):
                adapter.adapt(source, manifest, self.root / "out")
        self.assertFalse((self.root / "out").exists())

    def test_repeat_run_never_overwrites_previous_bytes(self):
        source, manifest = self.case()
        receipt = adapter.adapt(source, manifest, self.root / "out")
        before = {p.name: p.read_bytes() for p in (self.root / "out").iterdir()}
        with self.assertRaisesRegex(adapter.CompatibilityError, "already exists"):
            adapter.adapt(source, manifest, self.root / "out")
        self.assertEqual(before, {p.name: p.read_bytes() for p in (self.root / "out").iterdir()})

    def test_incomplete_source_run_rejected(self):
        source, manifest = self.case()
        (self.root / "INCOMPLETE.json").write_text("{}")
        with self.assertRaisesRegex(adapter.CompatibilityError, "INCOMPLETE"):
            adapter.adapt(source, manifest, self.root / "out")

    def test_both_real_r005_outputs_keep_561_endpoint_pairs_and_all_cells(self):
        if not os.environ.get('LAGCHECK_R005_CANONICAL_QA'):
            self.skipTest('Private R005 fixture: set LAGCHECK_R005_CANONICAL_QA to the archived canonical _qa directory')
        for mode in ("real_rate_only", "real_affine"):
            with self.subTest(mode=mode):
                origin = Path(os.environ['LAGCHECK_R005_CANONICAL_QA']) / mode
                source = origin / "corrected_laser_log.csv"
                before = source.read_bytes()
                result = adapter.adapt(source, origin / "correction_manifest.json", self.root / mode)
                self.assertEqual(result["firing_interval_count"], 561)
                self.assertEqual(result["output_row_count"], 3368)
                self.assertEqual(result["selected_sample"], "S10618 3 Down")
                self.assertEqual(result["metadata_records_with_changed_inherited_speed"], [2])
                self.assertEqual(Decimal(result["max_timestamp_quantization_error_s"]), Decimal("0.000500"))
                self.assertTrue(result["firing_endpoint_source_pairs_equal"])
                self.assertEqual(source.read_bytes(), before)
                output = Path(result["output_path"])
                rows = read_csv(output)[1:]
                original = read_csv(source)[1:]
                indices = result["source_record_by_output_record"]
                self.assertEqual(sorted(indices), list(range(1, 3369)))
                for new, old_id in zip(rows, indices):
                    self.assertEqual(new[1:], original[old_id - 1][1:])
                self.assertTrue(all(a[0] <= b[0] for a, b in zip(rows, rows[1:])))
                manifest = json.loads((output.parent / "correction_manifest.json").read_text(encoding="utf-8"))
                self.assertEqual(manifest["output_sha256"], hashlib.sha256(output.read_bytes()).hexdigest())
                self.assertFalse(manifest["full_export_diagnostic_passed"])
                self.assertFalse(manifest["importer_roundtrip_tested"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
