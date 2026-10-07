"""Semantic regressions for documented zero events and inherited membership."""
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

from lagcheck import iolite_adapter as adapter
from lagcheck.laser import load_laser, export_corrected


HEADER = ['Timestamp', 'Comment', 'Sequence Number', 'Scan Velocity (um/s)', 'Laser State', 'X']


def payload(rows):
    stream = io.StringIO(newline='')
    writer = csv.writer(stream)
    writer.writerow(HEADER)
    writer.writerows(rows)
    return stream.getvalue().encode('utf8')


def row(second, state, sample='', sequence='', speed='', x='0'):
    return ['2026-01-01 00:00:' + second, sample, sequence, speed, state, x]


class AdapterR23Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=Path(__file__).parent)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def fixture(self):
        original = self.root / 'raw.csv'
        original.write_bytes(payload([row('00.000000', 'Off', 'A', '1', '10'),
            row('01.000000', 'On'), row('02.000000', 'Off'),
            row('03.000000', 'On'), row('03.000000', 'Off'),
            row('04.000000', 'On'), row('05.000000', 'Off')]))
        options = dict(format='laser_state_csv', time_kind='datetime', sample_selector='A',
                       expected_sequence='1', zero_duration_policy='exclude')
        laser = load_laser(original, options)
        canonical = self.root / 'provisional_laser_log.csv'
        receipt = export_corrected(laser, canonical, 170.)
        manifest = dict(provisional=True, verified=False, mode='rate_only', rate_ppm=170.,
                        full_export_diagnostic_passed=False, output_sha256=adapter.sha256(canonical.read_bytes()),
                        input_sha256=dict(laser=adapter.sha256(original.read_bytes()), ms='not-read-by-adapter'),
                        export_plan=dict(mode='rate_only', offset_seconds='0'))
        report = dict(input_metadata=dict(laser=laser.metadata), input_sha256=manifest['input_sha256'],
                      provisional_export_verification=receipt)
        self.save('report.json', report)
        self.save('config.json', dict(laser=dict(options, path=str(original))))
        manifest_path = self.save('provisional_manifest.json', manifest)
        return original, canonical, manifest_path

    def save(self, name, data):
        path = self.root / name
        path.write_text(json.dumps(data), encoding='utf8')
        return path

    def test_documented_zero_is_omitted_with_exact_trace_and_positive_pairs_preserved(self):
        original, canonical, manifest = self.fixture()
        before = {p: p.read_bytes() for p in (original, canonical, manifest)}
        result = adapter.adapt(canonical, manifest, self.root / 'out', allow_provisional=True)
        self.assertEqual((result['rows_removed'], result['output_row_count'], result['firing_interval_count']), (2, 5, 2))
        self.assertEqual(result['firing_endpoint_source_records'], [[2, 3], [6, 7]])
        audit = json.loads((self.root / 'out/excluded_zero_duration_records.json').read_text())
        canonical_rows = list(csv.reader(io.StringIO(canonical.read_text())))[1:]
        self.assertEqual(audit['records'], [dict(source_record=i, cells=canonical_rows[i-1]) for i in (4, 5)])
        self.assertEqual(result['source_record_by_output_record'], [1, 2, 3, 6, 7])
        self.assertEqual(before, {p: p.read_bytes() for p in before})
        self.assertFalse(result['importer_roundtrip_tested'])

    def test_zero_without_recorded_report_is_rejected(self):
        _, canonical, manifest = self.fixture()
        self.save('report.json', {})
        with self.assertRaisesRegex(adapter.CompatibilityError, 'nonpositive'):
            adapter.adapt(canonical, manifest, self.root / 'out', allow_provisional=True)

    def test_forged_excluded_original_pair_is_rejected(self):
        _, canonical, manifest = self.fixture()
        path = self.root / 'report.json'
        report = json.loads(path.read_text())
        report['input_metadata']['laser']['input_processing']['excluded_intervals'][0]['start_record'] = 2
        self.save('report.json', report)
        with self.assertRaisesRegex(adapter.CompatibilityError, 'disagrees with original input'):
            adapter.adapt(canonical, manifest, self.root / 'out', allow_provisional=True)

    def test_changed_raw_bytes_are_rejected(self):
        original, canonical, manifest = self.fixture()
        original.write_bytes(original.read_bytes() + b'\n')
        with self.assertRaisesRegex(adapter.CompatibilityError, 'Original laser input changed'):
            adapter.adapt(canonical, manifest, self.root / 'out', allow_provisional=True)

    def test_reject_policy_cannot_authorize_zero_removal(self):
        _, canonical, manifest = self.fixture()
        config = json.loads((self.root / 'config.json').read_text())
        config['laser']['zero_duration_policy'] = 'reject'
        self.save('config.json', config)
        with self.assertRaisesRegex(adapter.CompatibilityError, 'policy/count'):
            adapter.adapt(canonical, manifest, self.root / 'out', allow_provisional=True)

    def test_sort_crossing_first_label_materializes_only_membership(self):
        data = payload([row('01.000000', 'Off', ' A ', '1', '10', ' 01 '),
                        row('00.000000', 'Off', x='retain'),
                        row('02.000000', 'On'), row('03.000000', 'Off')])
        output, mapping, details = adapter.prepare(data)
        rows = list(csv.reader(io.StringIO(output.decode())))[1:]
        self.assertEqual(rows[0][1:3], ['A', '1'])
        self.assertEqual(rows[0][-1], 'retain')
        self.assertEqual(rows[1][1], ' A ')
        self.assertEqual([m['source_record'] for m in mapping], [2, 1, 3, 4])
        self.assertTrue(details['sample_and_sequence_inheritance_equal_all_records'])
        self.assertTrue(details['non_membership_cells_equal_by_source_record'])

    def test_distinct_sample_or_sequence_is_not_hidden_by_materialization(self):
        for sample, sequence in [('B', ''), ('', '2')]:
            with self.subTest(sample=sample, sequence=sequence):
                with self.assertRaisesRegex(adapter.CompatibilityError, 'Only one selected'):
                    adapter.prepare(payload([row('01.000000', 'Off', 'A', '1', '10'),
                        row('00.000000', 'Off', sample, sequence),
                        row('02.000000', 'On'), row('03.000000', 'Off')]))

    def test_unrecorded_zero_negative_and_precision_collapse_are_rejected(self):
        for start, end in [('01.000000', '01.000000'), ('02.000000', '01.000000'),
                           ('01.000100', '01.000400')]:
            with self.subTest(start=start, end=end):
                with self.assertRaisesRegex(adapter.CompatibilityError, 'nonpositive'):
                    adapter.prepare(payload([row('00.000000', 'Off', 'A', '1', '10'),
                                             row(start, 'On'), row(end, 'Off')]))

    def test_zero_authorization_cannot_remove_positive_pair(self):
        with self.assertRaisesRegex(adapter.CompatibilityError, 'do not match actual firing pairs'):
            adapter.prepare(payload([row('00.000000', 'Off', 'A', '1', '10'),
                                     row('01.000000', 'On'), row('02.000000', 'Off')]),
                            [dict(start_record=2, end_record=3)])

    def test_invalid_calendar_date_and_missing_firing_off_are_rejected(self):
        with self.assertRaisesRegex(adapter.CompatibilityError, 'Invalid timestamp'):
            adapter.timestamp_seconds('2026-02-30 00:00:00.000')
        with self.assertRaisesRegex(adapter.CompatibilityError, 'complete nonempty'):
            adapter.prepare(payload([row('00.000000', 'Off', 'A', '1', '10'), row('01.000000', 'On')]))


if __name__ == '__main__':
    unittest.main()
