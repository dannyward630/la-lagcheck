"""Lines at the ends of the signal record are left out of the blocks by default."""
import unittest
import numpy as np
from lagcheck.estimator import diagnose, record_edge_lines
from lagcheck.models import SignalData, LaserData, InputError
from test_estimator import synthetic
from test_log_screen_and_options import record_starting_at_first_signal


class RecordEdgeLineCountTests(unittest.TestCase):
    def lines(self, starts, length=10.):
        return [{'Start Time': float(s), 'End Time': float(s) + length} for s in starts]

    def test_counts_lines_within_the_guard_of_either_end(self):
        t = np.arange(0, 10000, .1)
        starts = [0.1, 2.9, 3.1, 30.] + list(np.arange(60., 9900., 20.)) + [9960., 9988.]
        edge = record_edge_lines(t, self.lines(starts), 3., 10)
        self.assertEqual(edge['guard_s'], 3.)
        # 0.1 and 2.9 s start within 3 s of the first sample; only the last line ends within 3 s of the last sample.
        self.assertEqual((edge['start'], edge['end'], edge['applied']), (2, 1, True))

    def test_guard_follows_the_setting(self):
        t = np.arange(0, 10000, .1)
        starts = [0.1, 4., 45.] + list(np.arange(60., 9900., 20.))
        self.assertEqual(record_edge_lines(t, self.lines(starts), 3., 10)['start'], 1)
        self.assertEqual(record_edge_lines(t, self.lines(starts), 5., 10)['start'], 2)

    def test_not_applied_when_too_few_lines_would_remain(self):
        t = np.arange(0, 3000, .1)
        starts = [0.1] + list(np.arange(40., 40. + 20. * 119, 20.))
        self.assertEqual(len(starts), 120)
        edge = record_edge_lines(t, self.lines(starts), 3., 10)
        self.assertEqual((edge['start'], edge['end'], edge['applied'], edge['available_start']), (0, 0, False, 1))

    def test_nothing_to_leave_out_when_the_record_has_margins(self):
        t = np.arange(0, 3000, .1)
        edge = record_edge_lines(t, self.lines(np.arange(22., 2800., 20.)), 3., 10)
        self.assertEqual((edge['start'], edge['end'], edge['applied']), (0, 0, False))


class RecordEdgeDiagnosisTests(unittest.TestCase):
    def test_record_starting_at_the_first_firing_is_supported_by_default(self):
        ms, laser = record_starting_at_first_signal()
        r, blocks = diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': 3.})
        self.assertEqual(r['status'], 'supported', r.get('failed_checks'))
        self.assertTrue(r['settings']['skip_record_edge_lines'])
        self.assertTrue(r['record_edge_lines']['enabled'])
        self.assertGreaterEqual(r['record_edge_lines']['skipped_at_start'], 1)
        self.assertEqual(r['record_boundary_blocks']['baseline'], [])
        self.assertEqual(r['record_boundary_blocks']['closure'], [])
        self.assertAlmostEqual(r['rate_ppm'], 100., delta=3.)
        # The first block is searched over the full range on both sides.
        self.assertTrue(blocks[0]['retained'])

    def test_switching_the_rule_off_restores_the_search_edge_failure(self):
        ms, laser = record_starting_at_first_signal()
        r, _ = diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': 3., 'skip_record_edge_lines': False})
        self.assertEqual(r['status'], 'qc_failed')
        self.assertIn('peak_not_at_search_boundary', r['failed_checks'])
        self.assertEqual((r['record_edge_lines']['enabled'], r['record_edge_lines']['skipped_at_start']), (False, 0))

    def test_record_with_margins_is_analysed_once_with_every_line(self):
        ms, laser = synthetic(100., 3.)
        r, _ = diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': 3.})
        off, _ = diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': 3., 'skip_record_edge_lines': False})
        self.assertEqual(r['status'], 'supported')
        self.assertFalse(r['record_edge_lines']['applied'])
        self.assertEqual((r['record_edge_lines']['skipped_at_start'], r['record_edge_lines']['skipped_at_end']), (0, 0))
        self.assertEqual(r['rate_ppm'], off['rate_ppm'])

    def test_lines_are_left_out_only_at_the_end_that_needs_it(self):
        ms, laser = record_starting_at_first_signal()
        r, _ = diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': 3.})
        self.assertTrue(r['record_edge_lines']['applied'])
        self.assertEqual(r['record_edge_lines']['sides'], ['start'])
        self.assertEqual(r['record_edge_lines']['skipped_at_end'], 0)
        self.assertEqual(r['record_edge_lines']['guard_s'], 30.)

    def test_setting_must_be_boolean(self):
        ms, laser = synthetic(100., 3.)
        with self.assertRaises(InputError):
            diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': 3., 'skip_record_edge_lines': 1})

    def test_coverage_is_measured_against_all_lines(self):
        ms, laser = record_starting_at_first_signal()
        r, _ = diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': 3.})
        off, _ = diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': 3., 'skip_record_edge_lines': False})
        self.assertEqual(r['duration_s'], off['duration_s'])
        self.assertLess(r['coverage_fraction'], 1.)
        self.assertGreater(r['coverage_fraction'], .95)

    def test_offset_is_not_rate_with_the_rule_on(self):
        rates = []
        for offset in (2., 8.):
            ms, laser = synthetic(100., offset)
            r, _ = diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': offset})
            self.assertEqual(r['status'], 'supported', r)
            rates.append(r['rate_ppm'])
        # Different lines can enter the blocks for the two offsets, so the agreement is looser than with the rule off.
        self.assertAlmostEqual(*rates, delta=1.5)

    def test_step_is_still_rejected(self):
        ms, laser = synthetic(step=2.)
        r, _ = diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': 3.})
        self.assertNotEqual(r['status'], 'supported')

    def test_log_running_past_the_record_cannot_claim_coverage(self):
        rng = np.random.default_rng(148)
        lines, cursor = [], 5.
        for _ in range(240):
            length = round(float(rng.uniform(5, 18)), 1)
            lines.append({'Start Time': cursor, 'End Time': round(cursor + length, 1)})
            cursor = round(cursor + length + float(rng.uniform(3, 10)), 1)
        lines[-1]['End Time'] += 3000.
        t = np.arange(0, cursor + 30, .1)
        signal = np.full(len(t), 2.)
        for row in lines:
            signal[(t >= round(row['Start Time'] - 3, 1)) & (t < round(row['End Time'] - 3, 1))] = 1000.
        visible = (t[-1] - t[0]) / (lines[-1]['End Time'] - lines[0]['Start Time'])
        r, _ = diagnose(SignalData(t, signal), LaserData(lines), {'tolerance_s': .2, 'alignment_offset_s': 0.})
        self.assertIn(r['status'], ('qc_failed', 'insufficient_evidence'))
        if 'coverage_fraction' in r:
            self.assertLessEqual(r['coverage_fraction'], visible + 1e-9)
            self.assertFalse(r['checks']['coverage'])

    def test_shifted_tail_still_contradicts_the_fit(self):
        rng = np.random.default_rng(148)
        lines, cursor = [], 5.
        for _ in range(240):
            length = round(float(rng.uniform(5, 18)), 1)
            lines.append({'Start Time': cursor, 'End Time': round(cursor + length, 1)})
            cursor = round(cursor + length + float(rng.uniform(3, 10)), 1)
        lines[-1]['End Time'] += 50.
        t = np.arange(0, cursor + 30, .1)
        signal = np.full(len(t), 2.)
        for i, row in enumerate(lines):
            lag = -11. if i >= 216 else -3.
            signal[(t >= round(row['Start Time'] + lag, 1)) & (t < round(row['End Time'] + lag, 1))] = 1000.
        r, _ = diagnose(SignalData(t, signal), LaserData(lines), {'tolerance_s': .2, 'alignment_offset_s': 0.})
        self.assertEqual(r['status'], 'qc_failed')
        self.assertTrue({'no_contradictory_high_quality_blocks', 'baseline_heldout', 'all_sensitivity_evidence_supported'} & set(r['failed_checks']))


if __name__ == '__main__':
    unittest.main()
