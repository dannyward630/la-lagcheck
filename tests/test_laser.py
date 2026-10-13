import csv
import io
import tempfile
import unittest
from pathlib import Path

from lagcheck.laser import load_laser, export_corrected
from lagcheck.models import InputError


class LaserTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def source(self, text, name="input.csv"):
        path = self.root / name
        path.write_bytes(text.encode("utf-8"))
        return path

    def numeric(self, text):
        return load_laser(self.source(text), {"time_kind": "numeric", "time_unit": "s"})

    def test_positive_and_negative_rate_keep_input_and_other_cells(self):
        text = "Timestamp, Laser State, Comment, X\n0,Off,A, 01 \n10,On,,4\n20,Off,,5\n30,On,,6\n40,Off,,7\n"
        source = self.source(text)
        before = source.read_bytes()
        for rate in (100000, -100000):
            data = load_laser(source, {"time_kind": "numeric", "time_unit": "s"})
            report = export_corrected(data, self.root / f"out{rate}.csv", rate)
            parsed = load_laser(report["output_path"], report["output_config"])
            factor = 1 + rate / 1e6
            self.assertAlmostEqual(parsed.lines[0]["End Time"] - parsed.lines[0]["Start Time"], 10 * factor)
            self.assertTrue(report["input_unchanged"])
            self.assertTrue(report["non_time_cells_verified"])
            self.assertFalse(report["alignment_offset_applied"])
            self.assertEqual(source.read_bytes(), before)
            self.assertEqual(parsed.table["rows"][0][3], " 01 ")

    def test_metadata_timestamp_disorder_preserved_and_corrected(self):
        data = self.numeric("Timestamp,Laser State,Comment,X\n0,Off,A,1\n3,Off,,2\n1,Off,,3\n2,On,,4\n4,Off,,5\n")
        self.assertEqual(data.metadata["nonmonotonic_metadata_timestamp_count"], 1)
        report = export_corrected(data, self.root / "out.csv", 100000)
        parsed = load_laser(report["output_path"], report["output_config"])
        self.assertEqual([r[3] for r in parsed.table["rows"]], ["1", "2", "3", "4", "5"])
        self.assertAlmostEqual(float(parsed.table["times"][0][0]), -.2)
        self.assertAlmostEqual(float(parsed.table["times"][1][0]), 3.1)
        self.assertEqual(len(report["warnings"]), 1)

    def test_datetime_preamble_header_precision_and_timezone_independence(self):
        source = self.source("ESL metadata\r\nTimestamp, Laser State, Comment\r\n2026-01-01 00:00:00.000,Off,A\r\n2026-01-01 00:00:01.000,On,\r\n2026-01-01 00:00:02.001,Off,\r\n")
        data = load_laser(source, {"header_row": 1})
        self.assertEqual(data.lines[0]["Start Time"], 1.)
        report = export_corrected(data, self.root / "out.csv", 170)
        output = Path(report["output_path"]).read_bytes()
        self.assertTrue(output.startswith(b"ESL metadata\r\nTimestamp, Laser State, Comment\r\n"))
        self.assertIn(b"00:00:02.001170", output)
        self.assertLessEqual(report["max_timestamp_quantization_error_s"], .0000005)

    def test_submicrosecond_input_and_timezone_preserved(self):
        source = self.source("Timestamp,Laser State,Comment\n2026-01-01T00:00:00.000000123+02:00,On,A\n2026-01-01T00:00:01.000000123+02:00,Off,\n")
        data = load_laser(source, {})
        self.assertEqual(data.lines[0]["End Time"], 1.)
        report = export_corrected(data, self.root / "out.csv", 100)
        self.assertIn("+02:00", Path(report["output_path"]).read_text())
        self.assertLessEqual(report["max_timestamp_quantization_error_s"], .0000005)

    def test_interval_csv_numeric_units(self):
        source = self.source("a;b;sample;value\n1000;2000;A;raw\n3000;5000;A;v\n")
        config = {"format": "interval_csv", "delimiter": ";", "time_kind": "numeric", "time_unit": "ms", "start_column": "a", "end_column": "b", "sample_column": "sample"}
        data = load_laser(source, config)
        self.assertEqual(data.lines[0]["Start Time"], 0.)
        self.assertEqual(data.lines[1]["End Time"], 4.)
        report = export_corrected(data, self.root / "out.csv", -500000)
        result = load_laser(report["output_path"], report["output_config"])
        self.assertEqual(result.lines[1]["End Time"], 2.)

    def test_numeric_time_requires_explicit_unit(self):
        source = self.source("Timestamp,Laser State,Comment\n0,On,A\n1,Off,\n")
        with self.assertRaises(InputError):
            load_laser(source, {"time_kind": "numeric"})

    def test_multisample_requires_selector_and_keeps_selected_metadata(self):
        source = self.source("Timestamp,Laser State,Comment\n0,Off,A\n1,On,\n2,Off,\n3,Off,B\n4,On,\n5,Off,\n6,Off,\n")
        with self.assertRaises(InputError):
            load_laser(source, {"time_kind": "numeric", "time_unit": "s"})
        config = {"time_kind": "numeric", "time_unit": "s", "sample_selector": "B"}
        data = load_laser(source, config)
        self.assertEqual(len(data.table["rows"]), 4)
        self.assertEqual(data.lines[0]["Start Time"], 1.)
        report = export_corrected(data, self.root / "out.csv", 170)
        self.assertEqual(report["selected_row_count"], 4)

    def test_repeated_disjoint_name_and_selection_inside_interval_rejected(self):
        for text in (
            "Timestamp,Laser State,Comment\n0,On,A\n1,Off,\n2,On,B\n3,Off,\n4,On,A\n5,Off,\n",
            "Timestamp,Laser State,Comment\n0,On,B\n1,Off,A\n2,On,\n3,Off,\n",
        ):
            with self.assertRaises(InputError):
                load_laser(self.source(text), {"time_kind": "numeric", "time_unit": "s", "sample_selector": "A"})

    def test_unclosed_double_on_overlap_and_bad_metadata_rejected(self):
        for body in ("0,On,A\n", "0,On,A\n1,On,\n2,Off,\n",
                     "0,On,A\n2,Off,\n1,On,\n3,Off,\n",
                     "0,On,A\n0,Off,\n", "bad,Off,A\n1,On,\n2,Off,\n"):
            with self.assertRaises(InputError):
                self.numeric("Timestamp,Laser State,Comment\n" + body)

    def test_rate_bounds_and_nonoverwrite(self):
        data = self.numeric("Timestamp,Laser State,Comment\n0,On,A\n1,Off,\n")
        for rate in (float("nan"), float("inf"), -1000000, -1000001):
            output = self.root / "invalid.csv"
            with self.assertRaises(InputError):
                export_corrected(data, output, rate)
            self.assertFalse(output.exists())
        target = self.source("KEEP", "out.csv")
        with self.assertRaises(InputError):
            export_corrected(data, target, 170)
        self.assertEqual(target.read_text(), "KEEP")

    def test_changed_input_rejected(self):
        data = self.numeric("Timestamp,Laser State,Comment\n0,On,A\n1,Off,\n")
        Path(data.metadata["source_path"]).write_text("changed")
        output = self.root / "out.csv"
        with self.assertRaises(InputError):
            export_corrected(data, output, 170)
        self.assertFalse(output.exists())

    def test_custom_datetime_format_gets_microseconds(self):
        source = self.source("Timestamp,Laser State,Comment\n01/02/2026 03:04:05,On,A\n01/02/2026 03:04:06,Off,\n")
        data = load_laser(source, {"datetime_format": "%d/%m/%Y %H:%M:%S"})
        report = export_corrected(data, self.root / "out.csv", 170)
        self.assertEqual(report["output_config"]["datetime_format"], "%d/%m/%Y %H:%M:%S.%f")
        self.assertIn("06.000170", Path(report["output_path"]).read_text())

    def test_datetime_overflow_rejects_without_output(self):
        source = self.source("Timestamp,Laser State,Comment\n9999-12-31 23:59:58,On,A\n9999-12-31 23:59:59,Off,\n")
        data = load_laser(source, {})
        output = self.root / "out.csv"
        with self.assertRaises(InputError):
            export_corrected(data, output, 2000000)
        self.assertFalse(output.exists())

    def test_precision_collapse_rejects_without_output(self):
        data = self.numeric("Timestamp,Laser State,Comment\n0,On,A\n0.000001,Off,\n")
        output = self.root / "out.csv"
        with self.assertRaises(InputError):
            export_corrected(data, output, -999999.999999)
        self.assertFalse(output.exists())

    def test_affine_offset_uses_seconds_even_for_millisecond_cells(self):
        source=self.source('start,end,value\n1000,2000,A\n3000,5000,B\n')
        config=dict(format='interval_csv',time_kind='numeric',time_unit='ms',start_column='start',end_column='end')
        data=load_laser(source,config)
        receipt=export_corrected(data,self.root/'out.csv',-500000,offset_seconds='2.5')
        out=load_laser(receipt['output_path'],receipt['output_config'])
        self.assertEqual([float(row[0]) for row in out.table['rows']],[3500.,4500.])
        self.assertEqual([float(row[1]) for row in out.table['rows']],[4000.,5500.])
        self.assertEqual([row[2] for row in out.table['rows']],['A','B'])
        self.assertTrue(receipt['alignment_offset_applied'])
        self.assertEqual(source.read_text(),'start,end,value\n1000,2000,A\n3000,5000,B\n')

    def test_affine_datetime_crosses_midnight_without_changing_timezone(self):
        source=self.source('Timestamp,Laser State,Comment\n2026-01-01T23:59:59+02:00,On,A\n2026-01-02T00:00:00+02:00,Off,\n')
        data=load_laser(source,{})
        receipt=export_corrected(data,self.root/'out.csv',170,offset_seconds=2)
        text=Path(receipt['output_path']).read_text()
        self.assertIn('2026-01-02T00:00:01.000000+02:00',text)
        self.assertIn('2026-01-02T00:00:02.000170+02:00',text)


if __name__ == "__main__":
    unittest.main()
