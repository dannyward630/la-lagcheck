import csv
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from lagcheck.cli import run as _run, screen_log, main
from lagcheck.estimator import diagnose, screen_line_changes
from lagcheck.models import SignalData, LaserData, InputError
from lagcheck.signal_selection import select_signal
from test_estimator import synthetic


def run(*args, **kwargs):
    kwargs.setdefault('preview', False)
    return _run(*args, **kwargs)


def regular_lines(count=120, length=40., gap=3., jitter=.05, seed=7):
    rng = np.random.default_rng(seed)
    lines, cursor = [], 5.
    for _ in range(count):
        lines.append({'Start Time': cursor, 'End Time': cursor + length})
        cursor += length + gap + float(rng.normal(0, jitter))
    return lines


def record_starting_at_first_signal(rate=100., offset=3.):
    """The signal record begins as the signal of the first firing arrives."""
    ms, laser = synthetic(rate, offset)
    keep = ms.time >= laser.lines[0]['Start Time'] + offset - 1e-9
    return SignalData(ms.time[keep], ms.signal[keep]), laser


class LineChangeScreenTests(unittest.TestCase):
    def test_regular_raster_is_not_flagged(self):
        result = screen_line_changes(regular_lines())
        self.assertTrue(result['evaluated'])
        self.assertEqual(result['flagged'], [])
        self.assertLess(result['largest_unflagged_abs_deviation_s'], 1.)

    def test_excess_and_short_intervals_are_located(self):
        lines = regular_lines()
        for row in lines[48:]:
            row['Start Time'] += 60.2; row['End Time'] += 60.2
        for row in lines[90:]:
            row['Start Time'] -= 2.; row['End Time'] -= 2.
        result = screen_line_changes(lines)
        self.assertEqual([f['after_interval'] for f in result['flagged']], [48, 90])
        self.assertAlmostEqual(result['flagged'][0]['deviation_s'], 60.2, delta=.3)
        self.assertAlmostEqual(result['flagged'][1]['deviation_s'], -2., delta=.3)

    def test_step_below_the_limit_is_not_claimed(self):
        lines = regular_lines(jitter=.01)
        for row in lines[48:]:
            row['Start Time'] += .6; row['End Time'] += .6
        self.assertEqual(screen_line_changes(lines)['flagged'], [])
        self.assertEqual([f['after_interval'] for f in screen_line_changes(lines, limit_s=.4)['flagged']], [48])

    def test_local_scatter_bounds_the_sensitivity(self):
        # With 0.05 s scatter the threshold is ten scaled median absolute
        # deviations (about 0.7 s), whatever smaller limit is requested.
        lines = regular_lines(jitter=.05)
        for row in lines[48:]:
            row['Start Time'] += .4; row['End Time'] += .4
        self.assertEqual(screen_line_changes(lines, limit_s=.1)['flagged'], [])

    def test_irregular_intervals_are_not_flagged_throughout(self):
        _, laser = synthetic()                       # gaps drawn from 3 to 10 s
        self.assertEqual(screen_line_changes(laser.lines)['flagged'], [])
        lines = [dict(row) for row in laser.lines]
        for row in lines[100:]:
            row['Start Time'] += 60.; row['End Time'] += 60.
        self.assertEqual([f['after_interval'] for f in screen_line_changes(lines)['flagged']], [100])

    def test_too_few_intervals_are_reported_as_not_evaluated(self):
        result = screen_line_changes(regular_lines(count=3))
        self.assertFalse(result['evaluated'])
        self.assertEqual(result['flagged'], [])


class HiddenStepTests(unittest.TestCase):
    def hidden_step(self):
        """A 60 s timestamp step inside the first block, as in a log whose early lines were stamped by a stepped clock."""
        ms, laser = synthetic(100., 3.)
        lines = [dict(row) for row in laser.lines]
        for row in lines[:20]:
            row['Start Time'] -= 60.; row['End Time'] -= 60.
        return ms, LaserData(lines)

    def test_six_checks_alone_do_not_see_a_step_inside_an_unusable_block(self):
        ms, laser = self.hidden_step()
        r, blocks = diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3., 'line_change_limit_s': 1e6})
        self.assertEqual(r['status'], 'supported', r.get('failed_checks'))
        self.assertFalse(blocks[0]['retained'])
        self.assertAlmostEqual(r['rate_ppm'], 100., delta=3.)

    def test_screen_withholds_supported_status_and_locates_the_step(self):
        ms, laser = self.hidden_step()
        r, _ = diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3.})
        self.assertEqual(r['status'], 'qc_failed')
        self.assertEqual(r['failed_checks'], ['no_flagged_line_change_interval'])
        self.assertEqual([f['after_interval'] for f in r['laser_log_screen']['flagged']], [20])
        self.assertAlmostEqual(r['laser_log_screen']['flagged'][0]['deviation_s'], 60., delta=4.)
        self.assertAlmostEqual(r['rate_ppm'], 100., delta=3.)      # the estimate stays available

    def test_screen_is_reported_when_a_prerequisite_stops_the_diagnosis(self):
        ms, laser = self.hidden_step()
        r, _ = diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3., 'n_blocks': 100})
        self.assertEqual(r['reason'], 'insufficient_laser_intervals')
        self.assertEqual(r['checks'], {})
        self.assertEqual(len(r['laser_log_screen']['flagged']), 1)

    def test_invalid_limit_is_rejected(self):
        ms, laser = synthetic(100., 3.)
        for value in (0, -1., True):
            with self.subTest(value=value), self.assertRaises(InputError):
                diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3., 'line_change_limit_s': value})


class RecordBoundaryOptionTests(unittest.TestCase):
    def test_default_keeps_the_search_edge_check(self):
        ms, laser = record_starting_at_first_signal()
        r, blocks = diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3.})
        self.assertEqual(r['status'], 'qc_failed')
        self.assertIn('peak_not_at_search_boundary', r['failed_checks'])
        self.assertEqual(r['record_boundary_blocks']['policy'], 'check')
        self.assertEqual([(b['block'], b['side']) for b in r['record_boundary_blocks']['baseline']], [(0, 'start')])
        self.assertEqual(r['record_boundary_blocks']['excluded_blocks'], [])
        self.assertTrue(blocks[0]['retained'])

    def test_declared_option_excludes_the_block_but_keeps_it_in_agreement(self):
        ms, laser = record_starting_at_first_signal()
        r, blocks = diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3., 'exclude_record_boundary_blocks': True})
        self.assertEqual(r['status'], 'supported_excluding_record_boundary_blocks', r.get('failed_checks'))
        self.assertEqual(r['failed_checks'], [])
        self.assertEqual(r['record_boundary_blocks']['excluded_blocks'], [0])
        self.assertNotIn(0, r['baseline']['retained'])
        self.assertIn('record_boundary_excluded_by_option', blocks[0]['reasons'])
        self.assertLess(abs(r['record_boundary_blocks']['baseline'][0]['lag_minus_fitted_line_s']), .2)
        self.assertAlmostEqual(r['rate_ppm'], 100., delta=3.)

    def test_excluded_block_that_disagrees_still_fails_agreement(self):
        ms, laser = record_starting_at_first_signal()
        # Delay the signal of the first block by 1 s: it still peaks at the
        # record start after alignment only if the option hid it, so shift the
        # log instead and keep the record start fixed.
        lines = [dict(row) for row in laser.lines]
        for row in lines[24:]:
            row['Start Time'] -= 1.; row['End Time'] -= 1.
        r, _ = diagnose(ms, LaserData(lines), {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3., 'exclude_record_boundary_blocks': True,
                                              'line_change_limit_s': 1e6})
        self.assertEqual(r['status'], 'qc_failed')
        self.assertTrue({'no_contradictory_high_quality_blocks', 'baseline_heldout'} & set(r['failed_checks']), r['failed_checks'])

    def test_option_must_be_boolean(self):
        ms, laser = synthetic(100., 3.)
        with self.assertRaises(InputError):
            diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3., 'exclude_record_boundary_blocks': 1})

    def test_option_changes_nothing_without_a_record_boundary_block(self):
        ms, laser = synthetic(100., 3.)
        base, _ = diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3.})
        opt, _ = diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3., 'exclude_record_boundary_blocks': True})
        self.assertEqual(base['status'], 'supported')
        self.assertEqual(opt['status'], 'supported')
        self.assertEqual(base['rate_ppm'], opt['rate_ppm'])


class RecordBoundaryWorkflowTests(unittest.TestCase):
    def write_case(self, root, **analysis):
        ms, laser = record_starting_at_first_signal()
        with (root/'signal.csv').open('w', newline='') as f:
            w = csv.writer(f); w.writerow(['time_ms', 'element'])
            w.writerows(zip(ms.time*1000, ms.signal))
        with (root/'laser.csv').open('w', newline='') as f:
            w = csv.writer(f); w.writerow(['start', 'end'])
            w.writerows((row['Start Time'], row['End Time']) for row in laser.lines)
        # Both readers rebase at their first value: the first firing then
        # coincides with the first stored sample, so the declared offset is zero.
        config = {'ms': {'path': 'signal.csv', 'format': 'csv', 'header_row': 0, 'encoding': 'utf-8', 'delimiter': ',',
                         'time_kind': 'numeric', 'time_unit': 'ms', 'time_column': 'time_ms', 'signal_columns': ['element']},
                  'laser': {'path': 'laser.csv', 'format': 'interval_csv', 'header_row': 0, 'encoding': 'utf-8',
                            'delimiter': ',', 'time_kind': 'numeric', 'time_unit': 's',
                            'start_column': 'start', 'end_column': 'end'},
                  'analysis': dict({'tolerance_s': .2, 'alignment_offset_s': 0., 'skip_record_edge_lines': False}, **analysis)}
        (root/'config.json').write_text(json.dumps(config), encoding='utf-8')

    def test_option_gives_a_qualified_status_and_only_a_provisional_log(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root = Path(folder); self.write_case(root, exclude_record_boundary_blocks=True)
            code, report = _run(root/'config.json', root/'out', True)      # default preview is drawn
            self.assertEqual((code, report['status']), (2, 'supported_excluding_record_boundary_blocks'))
            self.assertFalse(report['correction_exported'])
            self.assertFalse((root/'out/corrected_laser_log.csv').exists())
            self.assertTrue((root/'out/provisional_laser_log.csv').is_file())
            self.assertEqual(report['provisional_failed_checks'], [])
            self.assertIn('all timing checks passed, but only with the declared option', report['provisional_note'])
            self.assertNotIn('did not all pass', report['provisional_note'])
            self.assertEqual(report['provisional_readback']['status'], 'passed', report['provisional_readback'])
            self.assertEqual(report['preview']['status'], 'written')
            text = (root/'out/report.md').read_text(encoding='utf-8')
            self.assertIn('Blocks at a record boundary', text)
            self.assertIn('Policy: exclude', text)
            self.assertIn('Failed checks: none.', text)
            self.assertNotIn('the timing evidence did not pass', text)
            manifest = json.loads((root/'out/provisional_manifest.json').read_text(encoding='utf-8'))
            self.assertFalse(manifest['verified'])
            self.assertEqual(manifest['laser_log_screen']['flagged_intervals'], 0)
            self.assertEqual(manifest['laser_log_screen']['limit_s'], 1.0)

    def test_default_run_names_the_record_boundary_block(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root = Path(folder); self.write_case(root)
            code, report = run(root/'config.json', root/'out')
            self.assertEqual((code, report['status']), (2, 'qc_failed'))
            self.assertEqual(report['record_boundary_blocks']['policy'], 'check')
            self.assertIn('Policy: check', (root/'out/report.md').read_text(encoding='utf-8'))


class ReviewFindingsTests(unittest.TestCase):
    """Behaviour fixed or documented after the independent code review of 0.6.0."""

    def test_screen_rejects_invalid_arguments(self):
        for kwargs in ({'neighbours': 0}, {'neighbours': True}, {'limit_s': float('nan')}, {'scale_multiple': 0}):
            with self.subTest(kwargs=kwargs), self.assertRaises(InputError):
                screen_line_changes(regular_lines(), **kwargs)
        lines = regular_lines(); lines[5]['End Time'] = float('nan')
        with self.assertRaises(InputError):
            screen_line_changes(lines)

    def test_raised_limit_remains_visible(self):
        ms, laser = HiddenStepTests().hidden_step()
        r, _ = diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3., 'line_change_limit_s': 1e6})
        self.assertEqual(r['status'], 'supported')
        self.assertEqual(r['laser_log_screen']['flagged'], [])
        self.assertEqual(r['laser_log_screen']['n_flagged_at_default_limit'], 1)
        self.assertEqual(r['laser_log_screen']['default_limit_s'], 1.0)

    def test_block_without_a_usable_peak_is_not_counted_as_excluded(self):
        ms, laser = record_starting_at_first_signal()
        flat_until = laser.lines[30]['Start Time'] + 3. - 1.          # covers the first block of every block division
        signal = np.where(ms.time < flat_until, 2., ms.signal)
        base, _ = diagnose(SignalData(ms.time, signal), laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3.})
        option, blocks = diagnose(SignalData(ms.time, signal), laser,
                                  {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3., 'exclude_record_boundary_blocks': True})
        self.assertEqual(option['status'], base['status'])
        self.assertEqual(option['record_boundary_blocks']['excluded_blocks'], [])
        self.assertEqual(option['record_boundary_blocks']['excluded_in_sensitivity'], {})
        self.assertNotIn('record_boundary_excluded_by_option', blocks[0]['reasons'])

    def test_option_cannot_see_a_lag_beyond_the_record_boundary(self):
        # Documented limitation. The signal of the first block arrives 1 s before the line fitted by the other blocks,
        # partly before the record. Its measured lag sits at the boundary, so the agreement check does not see the
        # difference; only the lower correlation peak shows it. The default search-edge check rejects the same record.
        rng = np.random.default_rng(148)
        lines, cursor = [], 5.
        for _ in range(240):
            length = float(rng.uniform(5, 18))
            lines.append({'Start Time': cursor, 'End Time': cursor + length})
            cursor += length + float(rng.uniform(3, 10))
        t = np.arange(0, cursor + 30, .1)
        signal = np.ones(len(t)) * 2.
        origin = lines[0]['Start Time']
        for i, row in enumerate(lines):
            early = 1. if i < 24 else 0.
            start = origin + (1 + 100e-6) * (row['Start Time'] - origin) + 3. - early
            end = origin + (1 + 100e-6) * (row['End Time'] - origin) + 3. - early
            signal += np.maximum(0, np.minimum(t + .05, end) - np.maximum(t - .05, start)) / .1 * 1000
        signal += rng.normal(0, .05, len(t))
        keep = t >= lines[0]['Start Time'] + 3. - 1e-9                # the record starts at the line-consistent arrival
        ms, laser = SignalData(t[keep], signal[keep]), LaserData(lines)
        default, _ = diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3.})
        option, _ = diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3., 'exclude_record_boundary_blocks': True})
        self.assertIn('peak_not_at_search_boundary', default['failed_checks'])
        row = option['record_boundary_blocks']['baseline'][0]
        self.assertTrue(row['excluded_from_fit_by_option'])
        self.assertLess(abs(row['lag_minus_fitted_line_s']), .2)        # the 1 s difference is not visible in the lag
        self.assertTrue(option['checks']['no_contradictory_high_quality_blocks'])
        self.assertLess(row['peak'], row['median_peak_of_retained_blocks'] - .02)
        self.assertIn('beyond the record boundary', option['record_boundary_blocks']['limitation'])

    def test_stability_alternatives_are_listed_under_the_default_policy(self):
        ms, laser = record_starting_at_first_signal()
        r, _ = diagnose(ms, laser, {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3.})
        listing = r['record_boundary_blocks']
        self.assertEqual(listing['excluded_in_sensitivity'], {})
        self.assertIn('n_blocks=15', listing['in_sensitivity'])
        self.assertNotIn('all', listing['in_sensitivity'])

    def test_short_path_reports_the_screen_even_without_a_fit(self):
        lines = regular_lines(count=60)
        for row in lines[30:]:
            row['Start Time'] += 60.; row['End Time'] += 60.
        t = np.arange(0, lines[-1]['End Time'] + 30, .1)
        ms = SignalData(t, np.full(len(t), 2.) + np.random.default_rng(1).normal(0, .01, len(t)))
        r, _ = diagnose(ms, LaserData(lines), {'tolerance_s': .2, 'short_acquisition': True})
        self.assertNotIn('no_flagged_line_change_interval', r.get('checks', {}))
        state = r['short_acquisition']['check_execution']['no_flagged_line_change_interval']
        self.assertEqual((state['state'], state['source']), ('fail', 'laser_log_screen'))

    def test_selection_keeps_the_sum_on_a_tie(self):
        ms, laser = synthetic(100., 3.)
        one = SignalData(ms.time, ms.signal, {'signal_columns': ['a']})
        for scale in (1., 3.):
            with self.subTest(scale=scale):
                other = SignalData(ms.time, ms.signal * scale, {'signal_columns': ['b']})
                summed = SignalData(ms.time, ms.signal * (1 + scale), {'signal_columns': ['a', 'b']})
                chosen, record = select_signal(summed, laser, {'tolerance_s': .2}, lambda cols: {'a': one, 'b': other}[cols[0]])
                self.assertIs(chosen, summed)
                self.assertEqual(record['selected_columns'], ['a', 'b'])

    def test_selection_without_contrast_keeps_the_configured_signal(self):
        ms, laser = synthetic(100., 3.)
        flat = SignalData(ms.time, np.full(len(ms.time), 5.), {'signal_columns': ['a', 'b']})
        single = SignalData(ms.time, np.full(len(ms.time), 5.), {'signal_columns': ['a']})
        chosen, record = select_signal(flat, laser, {'tolerance_s': .2}, lambda cols: single)
        self.assertIs(chosen, flat)
        self.assertIn('No candidate showed temporal contrast', record['note'])
        with self.assertRaises(InputError):
            select_signal(flat, laser, {'tolerance_s': .2, 'coarse_window_s': float('inf')}, lambda cols: single)

    def test_screen_log_rejects_a_non_object_configuration_and_marks_nothing_evaluated(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root = Path(folder)
            (root/'bad.json').write_text('[]', encoding='utf-8')
            with self.assertRaises(InputError):
                screen_log(root/'bad.json')
            ScreenLogCommandTests().write_log(root, regular_lines(count=3))
            code, summary = screen_log(root/'config.json')
            self.assertEqual((code, summary['status'], summary['evaluated_segments']), (2, 'not_evaluated', 0))
            self.assertEqual(summary['unevaluated_segments'], len(summary['segments']))


class PackagedR01OptionTests(unittest.TestCase):
    def test_option_acting_only_in_an_alternative_setting_still_qualifies_the_status(self):
        # In the packaged R01 example the first block is a record-boundary block
        # only in the 12- and 15-block alternatives of the stability check.
        import importlib, sys
        folder = Path(__file__).resolve().parents[1]/'examples'/'R01'
        sys.path.insert(0, str(folder))
        try:
            r = importlib.import_module('reconstruction')
        finally:
            sys.path.remove(str(folder))
        keep = np.isfinite(r.y) & np.isfinite(r.t)
        ms = SignalData(r.t[keep], r.y[keep])
        preset, _ = diagnose(ms, LaserData(r.L), {'tolerance_s': .2, 'skip_record_edge_lines': False})
        option, _ = diagnose(ms, LaserData(r.L), {'tolerance_s': .2, 'skip_record_edge_lines': False, 'exclude_record_boundary_blocks': True})
        self.assertEqual(preset['failed_checks'], ['all_sensitivity_evidence_supported'])
        self.assertEqual(preset['laser_log_screen']['flagged'], [])
        self.assertEqual(option['status'], 'supported_excluding_record_boundary_blocks')
        self.assertEqual(option['record_boundary_blocks']['excluded_blocks'], [])
        self.assertEqual(option['record_boundary_blocks']['excluded_in_sensitivity'], {'n_blocks=12': [0], 'n_blocks=15': [0]})
        # The example's reference mapping is now the default diagnosis, so the rate of this option is stated here.
        self.assertAlmostEqual(option['rate_ppm'], 168.99125403067967, delta=1e-6)


class SignalSelectionTests(unittest.TestCase):
    def columns(self):
        ms, laser = synthetic(100., 3.)
        rng = np.random.default_rng(3)
        noisy = rng.lognormal(12., 1.5, len(ms.time))      # no firing contrast, dominates the sum
        clear = SignalData(ms.time, ms.signal, {'signal_columns': ['element']})
        noise = SignalData(ms.time, noisy, {'signal_columns': ['noise']})
        summed = SignalData(ms.time, ms.signal + noisy, {'signal_columns': ['element', 'noise']})
        return summed, clear, noise, laser

    def test_highest_alignment_peak_is_selected_and_declared(self):
        summed, clear, noise, laser = self.columns()
        chosen, record = select_signal(summed, laser, {'tolerance_s': .2},
                                       lambda cols: {'element': clear, 'noise': noise}[cols[0]])
        self.assertIs(chosen, clear)
        self.assertEqual(record['selected_columns'], ['element'])
        self.assertEqual([c['columns'] for c in record['candidates']], [['element', 'noise'], ['element'], ['noise']])
        peaks = {tuple(c['columns']): c['peak'] for c in record['candidates']}
        self.assertGreater(peaks[('element',)], peaks[('element', 'noise')])

    def test_single_column_keeps_the_configured_signal(self):
        _, clear, _, laser = self.columns()
        chosen, record = select_signal(clear, laser, {'tolerance_s': .2}, lambda cols: self.fail('no reload expected'))
        self.assertIs(chosen, clear)
        self.assertEqual(record['selected_columns'], ['element'])

    def write_case(self, root, selection=None):
        summed, clear, noise, laser = self.columns()
        with (root/'signal.csv').open('w', newline='') as f:
            w = csv.writer(f); w.writerow(['time_ms', 'element', 'noise'])
            w.writerows(zip(clear.time*1000, clear.signal, noise.signal))
        with (root/'laser.csv').open('w', newline='') as f:
            w = csv.writer(f); w.writerow(['start', 'end'])
            w.writerows((row['Start Time'], row['End Time']) for row in laser.lines)
        config = {'ms': {'path': 'signal.csv', 'format': 'csv', 'header_row': 0, 'encoding': 'utf-8', 'delimiter': ',',
                         'time_kind': 'numeric', 'time_unit': 'ms', 'time_column': 'time_ms',
                         'signal_columns': ['element', 'noise']},
                  'laser': {'path': 'laser.csv', 'format': 'interval_csv', 'header_row': 0, 'encoding': 'utf-8',
                            'delimiter': ',', 'time_kind': 'numeric', 'time_unit': 's',
                            'start_column': 'start', 'end_column': 'end'},
                  'analysis': {'tolerance_s': .2, 'alignment_offset_s': 8.}}
        if selection is not None:
            config['ms']['signal_selection'] = selection
        (root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        return config

    def test_cli_records_the_selection_and_saves_the_columns_used(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root = Path(folder); self.write_case(root, 'highest_alignment_peak')
            code, report = run(root/'config.json', root/'out')
            self.assertEqual(report['signal_selection']['selected_columns'], ['element'])
            saved = json.loads((root/'out/config.json').read_text(encoding='utf-8'))
            self.assertEqual(saved['ms']['signal_columns'], ['element'])
            # The saved config sums the columns actually used and does not select again.
            self.assertEqual(saved['ms']['signal_selection'], 'sum')
            self.assertEqual(report['input_metadata']['ms']['signal_columns'], ['element'])
            self.assertEqual((code, report['status']), (0, 'supported'))
            self.assertIn('Matching-signal selection', (root/'out/report.md').read_text(encoding='utf-8'))

    def test_cli_default_selects_by_alignment_peak(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root = Path(folder); self.write_case(root)
            _, report = run(root/'config.json', root/'out')
            self.assertEqual(report['signal_selection']['mode'], 'highest_alignment_peak')
            self.assertEqual(report['signal_selection']['selected_columns'], ['element'])
            self.assertEqual(report['input_metadata']['ms']['signal_columns'], ['element'])

    def test_cli_sum_keeps_the_configured_columns(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root = Path(folder); self.write_case(root, 'sum')
            _, report = run(root/'config.json', root/'out')
            self.assertNotIn('signal_selection', report)
            self.assertEqual(report['input_metadata']['ms']['signal_columns'], ['element', 'noise'])

    def test_saved_config_replays_the_same_matching_signal(self):
        # Rerunning the config written by a run must not select a different signal.
        for selection in ('sum', None, 'highest_alignment_peak'):
            with self.subTest(selection=selection), tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
                root = Path(folder); self.write_case(root, selection)
                _, first = run(root/'config.json', root/'out')
                saved = json.loads((root/'out/config.json').read_text(encoding='utf-8'))
                self.assertEqual(saved['ms']['signal_selection'], 'sum')
                self.assertEqual(saved['ms']['signal_columns'], first['input_metadata']['ms']['signal_columns'])
                (root/'replay.json').write_text(json.dumps(saved), encoding='utf-8')
                _, second = run(root/'replay.json', root/'out2')
                self.assertEqual(second['status'], first['status'])
                self.assertEqual(second.get('rate_ppm'), first.get('rate_ppm'))
                self.assertEqual(second['input_metadata']['ms']['signal_columns'], first['input_metadata']['ms']['signal_columns'])
                self.assertNotIn('signal_selection', second)

    def test_cli_rejects_an_unknown_selection_mode(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root = Path(folder); self.write_case(root, 'best')
            with self.assertRaises(InputError):
                run(root/'config.json', root/'out')
            self.assertFalse((root/'out').exists())


class ScreenLogCommandTests(unittest.TestCase):
    def write_log(self, root, lines):
        with (root/'laser.csv').open('w', newline='') as f:
            w = csv.writer(f); w.writerow(['start', 'end'])
            w.writerows((row['Start Time'], row['End Time']) for row in lines)
        (root/'config.json').write_text(json.dumps({'laser': {
            'path': 'laser.csv', 'format': 'interval_csv', 'header_row': 0, 'encoding': 'utf-8', 'delimiter': ',',
            'time_kind': 'numeric', 'time_unit': 's', 'start_column': 'start', 'end_column': 'end'}}), encoding='utf-8')

    def test_screen_log_needs_no_signal_file(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root = Path(folder); self.write_log(root, regular_lines())
            code, summary = screen_log(root/'config.json')
            self.assertEqual((code, summary['flagged_intervals']), (0, 0))
            self.assertEqual(summary['segments'][0]['n_firing_intervals'], 120)

    def test_screen_log_exit_code_marks_a_flagged_interval(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root = Path(folder); lines = regular_lines()
            for row in lines[48:]:
                row['Start Time'] += 60.2; row['End Time'] += 60.2
            self.write_log(root, lines)
            code, summary = screen_log(root/'config.json')
            self.assertEqual((code, summary['flagged_intervals']), (2, 1))
            self.assertEqual(summary['segments'][0]['screen']['flagged'][0]['after_interval'], 48)
            self.assertEqual(main(['screen-log', str(root/'config.json'), '--limit', '100']), 0)


if __name__ == '__main__':
    unittest.main()
