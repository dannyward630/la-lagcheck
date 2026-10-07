"""Signal adapters: synthetic edge cases and an attributed public LDR fixture."""

import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lagcheck.models import InputError
from lagcheck.ms import load_signal


class SignalImportTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent)
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "signal.csv"

    def write(self, content, encoding="utf-8"):
        self.path.write_bytes(content.encode(encoding))

    def generic(self, **changes):
        config = dict(format="csv", header_row=0, delimiter=",", encoding="utf-8",
                      time_column="time", time_kind="numeric", time_unit="s",
                      signal_columns=["Fe"])
        config.update(changes)
        return config

    def test_units_normalize_first_observation(self):
        for unit, values in (("s", [2, 3, 5]), ("ms", [2000, 3000, 5000]),
                             ("us", [2000000, 3000000, 5000000]),
                             ("min", [2 / 60, 3 / 60, 5 / 60])):
            with self.subTest(unit=unit):
                self.write("time,Fe\n" + "".join(f"{v},1\n" for v in values))
                data = load_signal(self.path, self.generic(time_unit=unit))
                np.testing.assert_allclose(data.time, [0, 1, 3])
                self.assertAlmostEqual(data.metadata["time_first_seconds"], 2)
                self.assertEqual(data.metadata["time_unit"], unit)

    def test_unknown_or_missing_numeric_units_rejected(self):
        self.write("time,Fe\n1,1\n2,2\n3,3\n")
        for unit in ("seconds", "ns", "", None, []):
            with self.subTest(unit=unit), self.assertRaises(InputError):
                load_signal(self.path, self.generic(time_unit=unit))
        config = self.generic()
        del config["time_unit"]
        with self.assertRaisesRegex(InputError, "time_unit"):
            load_signal(self.path, config)

    def test_explicit_mapping_and_metadata_not_summed(self):
        self.write("instrument header\nelapsed;Cr52;Fe56;temperature;operator\n"
                   "1000;2;3;990;A\n2000;4;5;991;B\n3000;6;7;992;C\n")
        data = load_signal(self.path, self.generic(
            header_row=1, delimiter=";", time_column="elapsed", time_unit="ms",
            signal_columns=["Cr52", "Fe56"]))
        np.testing.assert_array_equal(data.signal, [5, 9, 13])
        np.testing.assert_array_equal(data.time, [0, 1, 2])
        self.assertEqual(data.metadata["signal_columns"], ["Cr52", "Fe56"])
        self.assertEqual(data.metadata["row_count"], 3)
        self.assertEqual(data.metadata["time_first"], "1000")
        self.assertEqual(data.metadata["sha256"], hashlib.sha256(self.path.read_bytes()).hexdigest())
        self.assertEqual(data.metadata["source"], str(self.path.resolve()))

    def test_generic_requires_explicit_configuration(self):
        self.write("time,Fe\n1,1\n2,2\n3,3\n")
        for key in ("header_row", "delimiter", "encoding", "time_column",
                    "time_kind", "signal_columns", "format"):
            config = self.generic()
            del config[key]
            with self.subTest(key=key), self.assertRaises(InputError):
                load_signal(self.path, config)

    def test_generic_header_row_is_physical_after_multiline_preamble(self):
        for newline in ("\n", "\r\n"):
            with self.subTest(newline=repr(newline)):
                preamble = newline.join(['"instrument note', 'continued on next line"', ''])
                self.write(preamble + newline.join(["time,Fe", "10,2", "11,4", "13,8", ""]))
                data = load_signal(self.path, self.generic(header_row=2))
                np.testing.assert_array_equal(data.time, [0, 1, 3])
                np.testing.assert_array_equal(data.signal, [2, 4, 8])
                self.assertEqual(data.metadata["header_row"], 2)
                self.assertEqual(data.metadata["preamble_physical_text"], preamble)
                # A CSV record index of 1 is not the physical header location.
                with self.assertRaises(InputError):
                    load_signal(self.path, self.generic(header_row=1))

    def test_generic_preamble_is_not_parsed_and_errors_use_physical_rows(self):
        self.write('"unclosed non-CSV instrument note\nsecond preamble line\n'
                   'time,Fe\n10,2\n11,4\n13,8\n')
        data = load_signal(self.path, self.generic(header_row=2))
        np.testing.assert_array_equal(data.signal, [2, 4, 8])
        self.write('"instrument note\ncontinued"\ntime,Fe\n10,2\n11,invalid\n13,8\n')
        with self.assertRaisesRegex(InputError, "Row 5:"):
            load_signal(self.path, self.generic(header_row=2))
        with self.assertRaisesRegex(InputError, "header_row"):
            load_signal(self.path, self.generic(header_row=99))

    def test_nonfinite_missing_or_unsorted_data_rejected(self):
        cases = ["1,1\n1,2\n3,3\n", "1,1\n3,2\n2,3\n", "1,1\nnan,2\n3,3\n",
                 "1,1\ninf,2\n3,3\n", "1,1\n2,\n3,3\n", "1,1\n2,nan\n3,3\n",
                 "1,1\n2,inf\n3,3\n", "1,1\n,2\n3,3\n", "1,1\n\n2,2\n3,3\n"]
        for rows in cases:
            with self.subTest(rows=rows):
                self.write("time,Fe\n" + rows)
                with self.assertRaises(InputError):
                    load_signal(self.path, self.generic())

    def test_missing_or_duplicate_selected_channel_rejected(self):
        self.write("time,Fe\n1,1\n2,2\n3,3\n")
        for columns in (["Zn"], ["Fe", "Fe"], [], "Fe", ["time"], [None]):
            with self.subTest(columns=columns), self.assertRaises(InputError):
                load_signal(self.path, self.generic(signal_columns=columns))

    def test_duplicate_header_and_malformed_rows_rejected(self):
        for content in ("time,Fe,Fe\n1,1,1\n2,2,2\n3,3,3\n",
                        "time,Fe\n1,1\n2,2,extra\n3,3\n",
                        'time,Fe\n1,1\n2,"unfinished\n3,3\n'):
            with self.subTest(content=content):
                self.write(content)
                with self.assertRaises(InputError):
                    load_signal(self.path, self.generic())

    def test_too_few_rows_and_unknown_format_rejected(self):
        self.write("time,Fe\n1,1\n2,2\n")
        with self.assertRaisesRegex(InputError, "three"):
            load_signal(self.path, self.generic())
        with self.assertRaisesRegex(InputError, "Unsupported signal format"):
            load_signal(self.path, self.generic(format="some_vendor"))

    def test_datetime_iso_and_explicit_format(self):
        self.write("time,Fe\n2026-09-19T23:59:59.500,1\n"
                   "2026-09-20T00:00:00.000,2\n2026-09-20T00:00:01.500,3\n")
        config = self.generic(time_kind="datetime")
        del config["time_unit"]
        data = load_signal(self.path, config)
        np.testing.assert_array_equal(data.time, [0, .5, 2])
        self.assertIsNone(data.metadata["time_first_seconds"])
        self.write("time,Fe\n19/09/2026 12:00:00,1\n19/09/2026 12:00:01,2\n"
                   "19/09/2026 12:00:03,3\n")
        with self.assertRaises(InputError):
            load_signal(self.path, config)
        config["datetime_format"] = "%d/%m/%Y %H:%M:%S"
        np.testing.assert_array_equal(load_signal(self.path, config).time, [0, 1, 3])

    def test_datetime_timezone_offsets_and_mixed_awareness(self):
        self.write("time,Fe\n2026-09-19T12:00:00+01:00,1\n"
                   "2026-09-19T11:00:01+00:00,2\n2026-09-19T11:00:02Z,3\n")
        config = self.generic(time_kind="datetime")
        np.testing.assert_array_equal(load_signal(self.path, config).time, [0, 1, 2])
        self.write("time,Fe\n2026-09-19T12:00:00+01:00,1\n"
                   "2026-09-19T12:00:01,2\n2026-09-19T12:00:02,3\n")
        with self.assertRaisesRegex(InputError, "timezone"):
            load_signal(self.path, config)

    def test_nondefault_encoding(self):
        self.write("time,Fe,comment\n1,1,caf\u00e9\n2,2,caf\u00e9\n3,3,caf\u00e9\n", "cp1252")
        np.testing.assert_array_equal(load_signal(self.path, self.generic(encoding="cp1252")).signal,
                                      [1, 2, 3])
        with self.assertRaises(InputError):
            load_signal(self.path, self.generic())

    def test_masshunter_named_header_and_printed_footer(self):
        self.write("C:\\Example\\Batch.b\\001SMPL.d\nIntensity Vs Time,CPS\n"
                   "Acquired : example\nTime [Sec],C13,Fe56,temperature\n"
                   "0.0932,1,5,999\n0.1622,2,6,998\n0.2312,3,7,997\n"
                   "\n\n          Printed:06-Jul-26 14:11:48\n")
        data = load_signal(self.path, {"format": "agilent_masshunter_csv",
                                       "signal_columns": ["C13", "Fe56"]})
        np.testing.assert_allclose(data.time, [0, .069, .138])
        np.testing.assert_array_equal(data.signal, [6, 8, 10])
        self.assertEqual(data.metadata["header_row"], 3)
        self.assertEqual(data.metadata["trailing_blank_records"], 2)
        self.assertEqual(data.metadata["recognized_footer"], ["Printed:06-Jul-26 14:11:48"])
        with self.assertRaisesRegex(InputError, "isotope"):
            load_signal(self.path, {"format": "agilent_masshunter_csv",
                                    "signal_columns": ["temperature"]})

    def test_masshunter_without_marker_or_with_bad_rows_rejected(self):
        config = {"format": "agilent_masshunter_csv", "signal_columns": ["Fe56"]}
        for content in ("Time [Sec],Fe56\n1,1\n2,2\n3,3\n",
                        "Intensity Vs Time,CPS\nTime [Sec],Fe56\n1,1\n\n2,2\n3,3\n",
                        "Intensity Vs Time,CPS\nTime [Sec],Fe56\n1,1\n2,2\n3,3\n,\n",
                        "Intensity Vs Time,CPS\nTime [Sec],Fe56\n1,1\n2,2\n3,3\n"
                        "Printed:example\n4,4\n"):
            self.write(content)
            with self.subTest(content=content), self.assertRaises(InputError):
                load_signal(self.path, config)

    def test_masshunter_cannot_override_seconds(self):
        self.write("Intensity Vs Time,CPS\nTime [Sec],Fe56\n1,1\n2,2\n3,3\n")
        with self.assertRaisesRegex(InputError, "time_unit"):
            load_signal(self.path, {"format": "agilent_masshunter_csv",
                                    "signal_columns": ["Fe56"], "time_unit": "ms"})


class QtegraLDRTests(unittest.TestCase):
    fixture = Path(__file__).resolve().parent / "fixtures" / "thermo_qtegra_ldr_5points.csv"

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent)
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "modified_ldr.csv"
        self.content = self.fixture.read_text(encoding="utf-8")

    def config(self, **changes):
        config = {"format": "thermo_qtegra_ldr_csv", "time_unit": "s",
                  "signal_columns": ["153Eu", "182W"]}
        config.update(changes)
        return config

    def reject_modified(self, old, new):
        self.assertIn(old, self.content)
        self.path.write_text(self.content.replace(old, new, 1), encoding="utf-8")
        with self.assertRaises(InputError):
            load_signal(self.path, self.config())

    def test_public_five_point_fixture_and_provenance(self):
        expected_digest = "4878036640e2e15d537e2b9e73b550e161ec40a6d31a3b8df2cfd50a3f7da722"
        self.assertEqual(hashlib.sha256(self.fixture.read_bytes()).hexdigest(), expected_digest)
        data = load_signal(self.fixture, self.config())
        np.testing.assert_allclose(data.time, [0, 1.00496, 2.00990, 3.01483, 4.01978], atol=1e-12)
        np.testing.assert_allclose(data.signal, [0, 2.50000025000002, 2.50000025000002,
                                                10.0000020000004, 2.50000025000002])
        self.assertEqual(data.metadata["sha256"], expected_digest)
        self.assertEqual(data.metadata["row_count"], 5)
        self.assertEqual(data.metadata["header_row"], 13)
        self.assertEqual(data.metadata["time_first"], "0.20229")
        self.assertEqual(data.metadata["time_first_seconds"], .20229)
        self.assertEqual(data.metadata["time_unit_source"], "explicit_config")
        self.assertEqual(data.metadata["qtegra"]["software_version"], "2.7.2425.65")
        self.assertEqual(data.metadata["qtegra"]["machine"], "iCAP RQ")
        self.assertEqual(data.metadata["qtegra"]["channel_metadata"]["31P"]["dwell_time_raw"], "0.2")
        self.assertIn("not long-run", data.metadata["adapter_validation_scope"])
        self.assertEqual(data.metadata["trailing_blank_records"], 1)

    def test_unit_required_and_explicit_mapping_constraints(self):
        config = self.config()
        del config["time_unit"]
        with self.assertRaisesRegex(InputError, "time_unit"):
            load_signal(self.fixture, config)
        for changes in ({"time_unit": "seconds"}, {"delimiter": ";"},
                        {"time_kind": "datetime"}, {"time_column": "31P"},
                        {"header_row": 14}, {"signal_columns": ["Time"]},
                        {"signal_columns": ["Fe56"]}):
            with self.subTest(changes=changes), self.assertRaises(InputError):
                load_signal(self.fixture, self.config(**changes))
        data = load_signal(self.fixture, self.config(time_unit="ms", signal_columns=["31P"]))
        np.testing.assert_allclose(data.time, np.array([0, 1.00496, 2.00990, 3.01483, 4.01978]) / 1000)
        self.assertAlmostEqual(data.signal[0], 2010.16161699401)

    def test_unrecognized_vendor_version_machine_and_metadata_rejected(self):
        for old, new in (("Name=Qtegra", "Name=NotQtegra"),
                         ("Version=2.7.2425.65", "Version=unknown"),
                         ("File Version=1", "File Version=2"),
                         ("Machine=iCAP RQ", "Machine=iCAP TQ"),
                         ("11-18-2020", "13-18-2020"),
                         ("RF Generator:", "Unrecognized:"),
                         ("\n\nTime,", "\nTime,")):
            with self.subTest(old=old, new=new):
                self.reject_modified(old, new)

    def test_header_requires_single_time_and_trailing_empty_column(self):
        for old, new in (("Time,31P,153Eu,182W,", "Time,31P,153Eu,182W"),
                         ("Time,31P,153Eu,182W,", "Time,31P,153Eu,182W,extra"),
                         ("Time,31P,153Eu,182W,", "Time,Time,153Eu,182W,"),
                         ("Time,31P,153Eu,182W,", "Time,temperature,153Eu,182W,"),
                         ("Time,31P,153Eu,182W,", "Elapsed,31P,153Eu,182W,")):
            with self.subTest(new=new):
                self.reject_modified(old, new)

    def test_dwell_row_is_required_validated_and_never_signal(self):
        for old, new in (("dwell time=0.2", "dwell time="),
                         ("dwell time=0.2", "dwell time=-1"),
                         ("dwell time=0.2", "dwell time=nan"),
                         ("xcal factor=66159.13097", "xcal factor=inf"),
                         ("xcal factor=66159.13097", "xcal factor=0"),
                         (",dwell time=0.2", "123,dwell time=0.2"),
                         ("xcal factor=98840.91543,", "xcal factor=98840.91543")):
            with self.subTest(old=old, new=new):
                self.reject_modified(old, new)
        metadata_row = self.content.splitlines()[14]
        self.reject_modified(metadata_row + "\n", "")

    def test_ldr_strict_numeric_time_signal_and_row_shape(self):
        for old, new in (("1.20725,1885.14213971733,2.50000025000002,0,",
                         "0.20229,1885.14213971733,2.50000025000002,0,"),
                        ("1.20725,", "nan,"),
                        ("1.20725,", "0.1,"),
                        ("1.20725,1885.14213971733,2.50000025000002,0,",
                         "1.20725,1885.14213971733,,0,"),
                        ("1.20725,1885.14213971733,2.50000025000002,0,",
                         "1.20725,1885.14213971733,inf,0,"),
                        ("1.20725,1885.14213971733,2.50000025000002,0,",
                         "1.20725,1885.14213971733,2.50000025000002,0,unexpected"),
                        ("\n1.20725,", "\n\n1.20725,"),
                        ("1.20725,1885.14213971733,2.50000025000002,0,", ",,,,")):
            with self.subTest(old=old, new=new):
                self.reject_modified(old, new)

    def test_ldr_is_not_generic_csv_or_masshunter(self):
        for config in ({"format": "csv", "header_row": 13, "delimiter": ",",
                        "encoding": "utf-8", "time_column": "Time", "time_kind": "numeric",
                        "time_unit": "s", "signal_columns": ["31P"]},
                       {"format": "agilent_masshunter_csv", "signal_columns": ["31P"]}):
            with self.subTest(format=config["format"]), self.assertRaises(InputError):
                load_signal(self.fixture, config)


if __name__ == "__main__":
    unittest.main()
