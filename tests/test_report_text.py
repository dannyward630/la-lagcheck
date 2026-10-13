"""Explanatory text of report.md, of the printed result and of the error messages.

These texts describe a completed diagnosis. They read numbers that the report
already holds and never change a status, a check, a rate or an exit code.
"""
import contextlib
import csv
import inspect
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import warnings

import numpy as np

from examples.make_iolite_demo import make_demo
from lagcheck import cli, estimator
from lagcheck.cli import run as _run, main
from lagcheck.estimator import diagnose
from lagcheck.laser import load_laser
from lagcheck.models import InputError, SignalData
from lagcheck.ms import load_signal
import test_cli
from test_estimator import synthetic


def run(*args, **kwargs):
    kwargs.setdefault('preview', False)
    return _run(*args, **kwargs)


def markdown(report, blocks):
    """report.md of a diagnosis made without files, as the command line would write it."""
    report = dict(report, correction_exported=False)
    cli._annotate(report, None, False, True)
    return report, cli._markdown(report, blocks).splitlines()


class Case(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def case(self, name='case', step=0., lines=None, **analysis):
        """The interval-CSV fixture of test_cli (240 lines, 100 µs/s), optionally cut to its first lines."""
        root = self.root/name
        root.mkdir()
        config = test_cli.WorkflowTests.make_case(None, root, step)
        config['analysis'].update({'skip_record_edge_lines': False}, **analysis)
        (root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        if lines is not None:
            with (root/'laser.csv').open(newline='') as f:
                rows = list(csv.reader(f))
            with (root/'laser.csv').open('w', newline='') as f:
                csv.writer(f).writerows(rows[:lines+1])
        return root

    def text(self, out):
        return (Path(out)/'report.md').read_text(encoding='utf-8').splitlines()

    def command(self, *args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = main([str(a) for a in args])
        return code, out.getvalue(), err.getvalue()


class RateTextTests(Case):
    def test_standard_error_is_printed_with_the_single_pass_rate(self):
        root = self.case()
        code, r = run(root/'config.json', root/'out')
        self.assertEqual((code, r['status'], r['closure']['converged']), (0, 'supported', True))
        text = self.text(root/'out')
        self.assertIn(f"Rate estimate (µs s⁻¹): {r['rate_ppm']}", text)
        self.assertIn('Kind of rate estimate: rate after iteration, reported without a standard error', text)
        self.assertIn(f"Single-pass rate (µs s⁻¹): {r['baseline']['rate_ppm']}; standard error of its slope (µs s⁻¹): "
                      f"{r['baseline']['slope_se_ppm']} ({r['baseline']['n_retained']} blocks)", text)
        code, printed, _ = self.command('run', root/'config.json', '--output', root/'cli', '--no-preview')
        summary = json.loads(printed)
        self.assertEqual((code, summary['rate_kind'], summary['rate_ppm']), (0, 'iterated', r['rate_ppm']))
        self.assertEqual(summary['single_pass_rate_ppm'], r['baseline']['rate_ppm'])
        self.assertEqual(summary['single_pass_rate_se_ppm'], r['baseline']['slope_se_ppm'])
        self.assertEqual(summary['notes'], [])

    def test_rate_after_iteration_is_never_paired_with_a_standard_error(self):
        # With this convergence limit the iteration ends at its second trial rate, which differs from the single-pass rate.
        ms, laser = synthetic(170., 3.)
        r, blocks = diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': 3., 'skip_record_edge_lines': False,
                                        'closure_tolerance_ppm': .6})
        self.assertTrue(r['closure']['converged'])
        self.assertGreater(abs(r['rate_ppm']-r['baseline']['rate_ppm']), .1)
        r, text = markdown(r, blocks)
        rate_line = next(line for line in text if line.startswith('Rate estimate'))
        self.assertEqual(rate_line, f"Rate estimate (µs s⁻¹): {r['rate_ppm']}")
        se = str(r['baseline']['slope_se_ppm'])
        self.assertEqual([line for line in text if se in line],
                         [f"Single-pass rate (µs s⁻¹): {r['baseline']['rate_ppm']}; standard error of its slope (µs s⁻¹): "
                          f"{se} ({r['baseline']['n_retained']} blocks)"])
        self.assertFalse(any('±' in line or str(r['closure']['residual_fit']['slope_se_ppm']) in line for line in text))
        summary = cli._rate_summary(r)
        self.assertEqual(summary, dict(rate_kind='iterated', single_pass_rate_ppm=r['baseline']['rate_ppm'],
                                       single_pass_rate_se_ppm=r['baseline']['slope_se_ppm'],
                                       single_pass_blocks=r['baseline']['n_retained']))

    def test_unconverged_iteration_states_the_reported_rate_and_the_last_trial_rate(self):
        root = self.case(closure_tolerance_ppm=1e-9, closure_iterations=2)
        code, r = run(root/'config.json', root/'out')
        last = r['closure']['history'][-1]
        self.assertEqual((code, r['closure']['converged'], len(r['closure']['history'])), (2, False, 2))
        self.assertEqual(r['rate_ppm'], r['baseline']['rate_ppm'])          # which rate is reported is unchanged
        note = next(n for n in r['notes'] if n.startswith('The iteration did not converge.'))
        self.assertIn(f"the trial rate was {last['trial_ppm']:.3f} µs s⁻¹", note)
        self.assertIn(f"The reported rate estimate is the single-pass rate ({r['baseline']['rate_ppm']:.3f} µs s⁻¹), "
                      'not the last trial rate.', note)
        text = self.text(root/'out')
        self.assertIn('Kind of rate estimate: single-pass rate; the iteration did not converge', text)
        self.assertIn('- ' + note, text)
        detail = text[text.index('- closure_converged: not met')+1]
        self.assertTrue(detail.startswith('  - After 2 iteration(s) the trial rate was'), detail)
        self.assertIn('convergence needs an absolute residual below 1e-09', detail)
        code, printed, _ = self.command('run', root/'config.json', '--output', root/'cli', '--no-preview')
        summary = json.loads(printed)
        self.assertEqual((code, summary['rate_kind']), (2, 'single_pass'))
        self.assertEqual(summary['rate_ppm'], summary['single_pass_rate_ppm'])
        self.assertIn(note, summary['notes'])

    def test_no_rate_gives_no_rate_lines(self):
        root = self.case(peak_min=.999)
        _, r = run(root/'config.json', root/'out')
        self.assertIsNone(r['rate_ppm'])
        self.assertFalse(any(line.startswith(('Kind of rate estimate', 'Single-pass rate')) for line in self.text(root/'out')))
        self.assertEqual(cli._rate_summary(r)['rate_kind'], None)


class CheckTextTests(Case):
    def test_each_check_that_is_not_met_gives_its_measured_value_and_limit(self):
        root = self.case(step=2.)
        code, r = run(root/'config.json', root/'out')
        self.assertEqual((code, r['status']), (2, 'qc_failed'))
        text = self.text(root/'out')
        for check in r['failed_checks']:
            position = text.index(f'- {check}: not met')
            self.assertTrue(text[position+1].startswith('  - '), check)
        for check, passed in r['checks'].items():
            if passed:
                self.assertFalse(text[text.index(f'- {check}: pass')+1].startswith('  - '), check)
        heldout = text[text.index('- baseline_heldout: not met')+1]
        self.assertEqual(heldout, '  - Largest held-out prediction error of the single-pass fit: '
                         f"{r['baseline']['heldout_max_s']:.4g} s; the tolerance is 0.2 s.")
        off_line = text[text.index('- no_contradictory_high_quality_blocks: not met')+1]
        self.assertTrue(off_line.startswith('  - Lag differs from the fitted line by more than the tolerance of 0.2 s: '))
        listed = [part.split()[1] for part in off_line.split(': ', 1)[1].split(', ')]
        self.assertEqual([int(b) for b in listed], r['contradictory_high_quality_blocks'])
        differences = [float(part.split('(')[1].split()[0]) for part in off_line.split(': ', 1)[1].split(', ')]
        self.assertTrue(all(abs(d) > .2 for d in differences), differences)
        envelope = text[text.index('- sensitivity_envelope: not met')+1]
        self.assertIn(f"this spread is {r['sensitivity_envelope_s']:.4g} s; the tolerance is 0.2 s.", envelope)

    def test_read_back_check_that_is_not_met_gives_its_measured_value_and_limit(self):
        import test_export_verification
        root = self.root/'affine'
        root.mkdir()
        test_export_verification.ExportVerificationTests.case(None, root)
        code, r = run(root/'config.json', root/'out', True, 'affine', '0', 'Known synthetic coordinate')
        post = r['post_export_verification']
        self.assertEqual((code, r['input_diagnostic_status'], r['status']), (2, 'supported', 'export_qc_failed'))
        self.assertIn('residual_rate_below_threshold', post['failed_checks'])
        text = self.text(root/'out')
        largest, limit = cli._against(max(abs(post['residual_baseline_rate_ppm']), abs(post['residual_rate_ppm'])),
                                      post['residual_rate_limit_ppm'])
        self.assertEqual(text[text.index('- residual_rate_below_threshold: not met')+1],
                         f'  - Larger of the two residual rates in absolute value: {largest} µs s⁻¹; it must be below {limit} µs s⁻¹.')
        self.assertEqual(limit, '0.2')
        self.assertEqual(cli._readback_detail('written_bytes_unchanged', post, r), [])
        self.assertEqual(cli._readback_detail('heldout', {}, r), [])           # an incomplete record gets no text

    def test_alternatives_without_a_rate_are_listed(self):
        root = self.case(lines=143)
        code, r = run(root/'config.json', root/'out')
        self.assertEqual((code, r['n_intervals']), (2, 143))
        self.assertFalse(r['checks']['all_sensitivity_configs_estimable'])
        text = self.text(root/'out')
        position = text.index('- all_sensitivity_configs_estimable: not met')
        self.assertEqual(text[position+1:position+3], ['  - n_blocks=12: needs 144 firing intervals; 143 are in the fits',
                                                      '  - n_blocks=15: needs 180 firing intervals; 143 are in the fits'])
        missing = [v['config'] for v in r['sensitivity'] if v['rate_ppm'] is None]
        self.assertEqual(missing, ['n_blocks=12', 'n_blocks=15'])

    def test_diagnosis_that_stops_before_the_checks_gives_the_reason_and_points_to_the_blocks(self):
        root = self.case(peak_min=.999)
        code, r = run(root/'config.json', root/'out')
        self.assertEqual((code, r['status'], r['reason'], r['checks']), (2, 'insufficient_evidence', 'fewer_than_four_usable_blocks', {}))
        text = self.text(root/'out')
        position = text.index('## Input diagnostic checks')
        self.assertEqual(text[position+2], 'No check was evaluated: the diagnosis stopped before the checks (fewer_than_four_usable_blocks).')
        self.assertEqual(text[position+3], f"{r['gate_passing_blocks']} of {r['computed_blocks']} blocks passed the correlation gates "
                         '(peak above 0.999, margin above 0.01). A fit needs at least four such blocks that also lie inside the signal record.')
        self.assertEqual(text[position+4], 'blocks.csv lists the lag, correlation peak, margin and reasons of every block.')
        self.assertTrue((root/'out/blocks.csv').is_file())

    def test_explanations_never_raise_on_an_incomplete_report(self):
        partial = dict(status='qc_failed', correction_exported=False, rate_ppm=1., settings={'tolerance_s': .2},
                       checks=dict.fromkeys(('coverage', 'closure_heldout', 'sensitivity_envelope', 'closure_converged', 'unknown_check'), False))
        text = cli._markdown(partial, None).splitlines()
        self.assertIn('- coverage: not met', text)
        self.assertIn('- unknown_check: not met', text)
        empty = {}
        cli._annotate(empty, 170.9, True, True)
        self.assertEqual(empty['notes'][-1], 'The reference rate (170.9 µs s⁻¹) was not applied to a log.')
        self.assertEqual(cli._rate_summary({'rate_ppm': 1., 'baseline': 'not a mapping'})['rate_kind'], None)
        self.assertEqual(cli._stopped_before_checks({'reason': 'fewer_than_four_usable_blocks'}),
                         ['No check was evaluated: the diagnosis stopped before the checks (fewer_than_four_usable_blocks).'])

    def test_texts_use_the_block_counts_and_alignment_criteria_of_the_estimator(self):
        ms, laser = synthetic(100., 3.)
        r, _ = diagnose(ms, laser, {'tolerance_s': .2, 'alignment_offset_s': 3., 'skip_record_edge_lines': False})
        counts = tuple(int(v['config'].split('=')[1]) for v in r['sensitivity'] if v['config'].startswith('n_blocks='))
        self.assertEqual(counts, cli.STABILITY_BLOCK_COUNTS)
        # The sentence in cli.ALIGNMENT_CRITERIA repeats the two limits written in estimator._coarse.
        self.assertIn('values[k] > .25 and margin > .01', inspect.getsource(estimator._coarse))
        self.assertIn('a peak above 0.25, a margin above 0.01', cli.ALIGNMENT_CRITERIA)

    def test_block_without_a_usable_peak_raises_no_floating_point_warning(self):
        ms, laser = synthetic(100., 3.)
        keep = ms.time >= laser.lines[0]['Start Time'] + 3. - 1e-9
        time, signal = ms.time[keep], ms.signal[keep]
        signal = np.where(time < laser.lines[30]['Start Time'] + 2., 2., signal)       # flat over the first block
        settings = {'tolerance_s': .2, 'skip_record_edge_lines': False, 'alignment_offset_s': 3.}
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            r, blocks = diagnose(SignalData(time, signal), laser, settings)
        # The case is present: in the finer block divisions the first block has no peak, and its margin stays undefined.
        undefined = [v['config'] for v in r['sensitivity'] for b in v['blocks'] if b['block'] == 0 and b['margin'] is None]
        self.assertIn('n_blocks=15', undefined)
        self.assertFalse(blocks[0]['retained'])


class NoteTests(Case):
    def test_short_record_names_the_line_count_and_the_checks_that_cannot_be_made(self):
        root = self.case(lines=143)
        code, r = run(root/'config.json', root/'out')
        self.assertEqual((code, r['status']), (2, 'qc_failed'))
        self.assertEqual(r['notes'], [
            'Record too short for the complete analysis: 143 firing intervals. '
            'The complete analysis needs 180 (15 blocks of at least 12). '
            'The stability alternative(s) n_blocks=12, n_blocks=15 cannot be estimated. '
            'The checks all_sensitivity_configs_estimable and all_sensitivity_evidence_supported need every alternative '
            'and are therefore reported as not met.'])
        self.assertIn('- ' + r['notes'][0], self.text(root/'out'))
        code, printed, _ = self.command('run', root/'config.json', '--output', root/'cli', '--no-preview')
        self.assertEqual((code, json.loads(printed)['notes']), (2, r['notes']))
        # The checks themselves are as before: the two that need every alternative are not met, the others are.
        self.assertEqual(r['failed_checks'], ['all_sensitivity_configs_estimable', 'all_sensitivity_evidence_supported'])

    def test_record_too_short_for_a_fit_says_that_no_check_was_evaluated(self):
        root = self.case(lines=100)
        code, r = run(root/'config.json', root/'out')
        self.assertEqual((code, r['status'], r['reason']), (2, 'insufficient_evidence', 'insufficient_laser_intervals'))
        self.assertEqual(r['notes'], ['Record too short: 100 firing intervals. A fit with 10 blocks needs at least 120 '
                                      '(12 per block), so no rate was estimated and no check was evaluated.'])
        text = self.text(root/'out')
        self.assertIn('100 firing intervals; 10 blocks of at least 12 need 120.', text)
        self.assertFalse(any('blocks.csv' in line for line in text))       # no block was measured, so no file to point to

    def test_complete_record_has_no_length_note_and_left_out_lines_are_counted(self):
        root = self.case()
        _, r = run(root/'config.json', root/'out')
        self.assertEqual(r['notes'], [])
        report = dict(n_intervals=182, settings=dict(n_blocks=10), checks=dict(all_sensitivity_configs_estimable=False),
                      record_edge_lines=dict(skipped_at_start=3, skipped_at_end=0))
        note = cli._record_length_note(report)
        self.assertIn('182 firing intervals, of which 179 are in the fits after the lines at the ends of the record were left out. '
                      'The complete analysis needs 180 (15 blocks of at least 12).', note)
        self.assertIn('The stability alternative(s) n_blocks=15 cannot be estimated.', note)
        self.assertIsNone(cli._record_length_note(dict(report, record_edge_lines=dict(skipped_at_start=2, skipped_at_end=0))))

    def test_reference_rate_that_is_not_used_is_stated_with_the_reason(self):
        root = self.case(step=2.)
        code, r = run(root/'config.json', root/'out', True, reference_rate_ppm=999.)
        self.assertEqual(code, 2)
        self.assertEqual(r['provisional_rate_source'], 'acquisition_estimate')            # the order of precedence is unchanged
        self.assertAlmostEqual(r['provisional_rate_ppm'], r['rate_ppm'])
        use = r['reference_rate']
        self.assertEqual((use['supplied_ppm'], use['used'], use['reason']), (999., False, 'acquisition_estimate_available'))
        self.assertEqual(use['note'], f"The reference rate (999 µs s⁻¹) was not used: a rate was estimated from this acquisition "
                         f"({r['rate_ppm']:.3f} µs s⁻¹), and that estimate takes precedence. "
                         'A reference rate is applied only when the acquisition gives no estimate or its estimate is not eligible for provisional export.')
        self.assertIn(use['note'], r['notes'])
        self.assertIn('- ' + use['note'], self.text(root/'out'))
        code, printed, _ = self.command('run', root/'config.json', '--output', root/'cli', '--export', '--no-preview',
                                        '--reference-rate-ppm', '170.9')
        self.assertEqual(code, 2)
        self.assertTrue(any(n.startswith('The reference rate (170.9 µs s⁻¹) was not used:') for n in json.loads(printed)['notes']))

    def test_reference_rate_use_is_stated_in_every_other_case(self):
        root = self.case(peak_min=.999)
        expected = [('used', dict(export=True), True, 'no_acquisition_estimate', 'applied to the provisional log'),
                    ('no_export', dict(export=False), False, 'export_not_requested', 'not applied to a log because --export was not requested'),
                    ('off', dict(export=True, provisional=False), False, 'provisional_output_disabled', 'provisional output is disabled')]
        for name, options, used, reason, words in expected:
            _, r = run(root/'config.json', root/name, reference_rate_ppm=170.9, **options)
            use = r['reference_rate']
            self.assertEqual((use['used'], use['reason']), (used, reason), name)
            self.assertIn(words, use['note'], name)
            self.assertEqual(r.get('provisional_rate_source'), 'user_reference_rate' if used else None, name)
        _, r = run(root/'config.json', root/'none', True)
        self.assertNotIn('reference_rate', r)
        self.assertEqual(r['notes'], [])

    def test_report_without_export_names_the_export_option(self):
        root = self.case()
        run(root/'config.json', root/'plain')
        hint = [line for line in self.text(root/'plain') if line.startswith('No corrected file requested.')]
        self.assertEqual(len(hint), 1)
        self.assertIn('run again with --export and a new output directory', hint[0])
        self.assertIn('provisional (not verified) log is written when a rate is available', hint[0])
        run(root/'config.json', root/'exported', True)
        self.assertFalse(any('No corrected file requested.' in line for line in self.text(root/'exported')))


class InputMessageTests(Case):
    def masshunter(self):
        path = self.root/'signal.csv'
        path.write_text('Synthetic record\nIntensity Vs Time,CPS\nTime [Sec],C13,Fe56\n0,1,2\n1,2,3\n2,3,4\n', encoding='utf-8')
        return path

    def test_missing_signal_column_lists_the_columns_of_the_file(self):
        path = self.masshunter()
        with self.assertRaises(InputError) as caught:
            load_signal(path, dict(format='agilent_masshunter_csv', signal_columns=['C13', 'Si29']))
        self.assertEqual(str(caught.exception),
                         "Columns not found in signal file: ['Si29']. Columns in the file: ['Time [Sec]', 'C13', 'Fe56'].")

    def test_wrong_channel_prints_the_input_error_only(self):
        laser = self.root/'laser.csv'
        laser.write_text('Timestamp,Laser State,Comment\n2026-01-01 00:00:00,On,Synthetic A\n2026-01-01 00:00:01,Off,\n', encoding='utf-8')
        config = dict(ms=dict(path=str(self.masshunter()), format='agilent_masshunter_csv', signal_columns=['C13', 'Si29']),
                      laser=dict(path=str(laser), format='laser_state_csv', time_kind='datetime'), analysis=dict(tolerance_s=.2))
        (self.root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        code, printed, errors = self.command('run', self.root/'config.json', '--output', self.root/'out')
        self.assertEqual((code, printed), (1, ''))
        self.assertEqual([json.loads(line)['status'] for line in errors.splitlines()], ['input_or_tool_error'])
        self.assertIn("Columns in the file: ['Time [Sec]', 'C13', 'Fe56']", errors)
        self.assertNotIn('Partial-signal inspection', errors)
        # The folder of the failed run is marked as before; the inspection view exists for numeric CSV mappings only.
        self.assertEqual(sorted(p.name for p in (self.root/'out').iterdir()), ['INCOMPLETE.json', 'config.json'])

    def test_partial_inspection_of_a_numeric_csv_mapping_is_still_attempted(self):
        root = self.case()
        config = json.loads((root/'config.json').read_text(encoding='utf-8'))
        config['ms']['signal_columns'] = ['element', 'absent']
        (root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        code, _, errors = self.command('run', root/'config.json', '--output', root/'out')
        self.assertEqual(code, 1)
        self.assertIn('Partial-signal inspection column must occur exactly once: absent', errors)
        self.assertIn("Columns not found in signal file: ['absent']. Columns in the file: ['time_ms', 'element', 'not_signal'].", errors)

    def log(self, names):
        path = self.root/'laser.csv'
        rows = ['Timestamp,Laser State,Comment']
        for i, name in enumerate(names):
            rows += [f'2026-01-01 00:{i:02d}:00,On,{name}', f'2026-01-01 00:{i:02d}:30,Off,']
        path.write_text('\n'.join(rows)+'\n', encoding='utf-8')
        return path

    def test_unknown_sample_selector_lists_the_sample_names(self):
        path = self.log(['Synthetic A', 'Synthetic B'])
        options = dict(format='laser_state_csv', time_kind='datetime')
        with self.assertRaises(InputError) as caught:
            load_laser(path, dict(options, sample_selector='Synthetic C'))
        self.assertEqual(str(caught.exception), "Selected sample 'Synthetic C' was not found. "
                         "Sample names in column 'Comment': 'Synthetic A', 'Synthetic B'")
        with self.assertRaises(InputError) as caught:
            load_laser(path, options)
        self.assertEqual(str(caught.exception), "Multiple samples: set sample_column and sample_selector. "
                         "Sample names in column 'Comment': 'Synthetic A', 'Synthetic B'")
        self.assertEqual(len(load_laser(path, dict(options, sample_selector='Synthetic B')).lines), 1)

    def test_long_list_of_sample_names_is_shortened(self):
        path = self.log([f'Synthetic {i:02d}' for i in range(25)])
        with self.assertRaises(InputError) as caught:
            load_laser(path, dict(format='laser_state_csv', time_kind='datetime', sample_selector='absent'))
        self.assertTrue(str(caught.exception).endswith("'Synthetic 19' and 5 more"), str(caught.exception))


class CommandLineTests(Case):
    def titles(self, *args):
        """Run a command and return the figure titles it drew."""
        from matplotlib.figure import Figure
        with patch.object(Figure, 'suptitle', autospec=True) as drawn:
            code, printed, errors = self.command(*args)
        self.assertEqual(errors, '')
        return code, [call.args[1] for call in drawn.call_args_list]

    def test_no_sample_label_switch_leaves_the_label_out_of_the_title(self):
        root = make_demo(self.root/'demo')
        code, titles = self.titles('run', root/'config.json', '--output', self.root/'default')
        self.assertEqual(code, 0)
        self.assertEqual(titles[0].splitlines()[0], 'Synthetic A; signal')            # the default title is unchanged
        self.assertNotIn('sample_label_in_title', json.loads((self.root/'default/preview_pair.json').read_text(encoding='utf-8')))
        code, titles = self.titles('run', root/'config.json', '--output', self.root/'plain', '--no-sample-label')
        self.assertEqual(code, 0)
        self.assertEqual(titles[0].splitlines()[0], 'signal')
        self.assertFalse(any('Synthetic A' in title for title in titles))
        self.assertIs(json.loads((self.root/'plain/preview_pair.json').read_text(encoding='utf-8'))['sample_label_in_title'], False)
        code, titles = self.titles('preview', self.root/'default', '--name', 'second', '--no-sample-label')
        self.assertEqual((code, titles[0].splitlines()[0]), (0, 'signal'))
        code, titles = self.titles('preview', self.root/'default', '--name', 'third')
        self.assertEqual((code, titles[0].splitlines()[0]), (0, 'Synthetic A; signal'))
        # The switch changes the picture only: both runs report the same result.
        reports = [json.loads((self.root/name/'report.json').read_text(encoding='utf-8')) for name in ('default', 'plain')]
        self.assertEqual(reports[0], reports[1])

    def test_no_sample_label_switch_reaches_separately_processed_fragments(self):
        root = make_demo(self.root/'demo')
        with (root/'laser.csv').open('a', encoding='utf-8', newline='') as f:
            w = csv.writer(f)
            w.writerow(['2026-01-01 05:00:00', 'Other', '2', '10', 'Off', 'meta'])
            w.writerow(['2026-01-01 06:00:00', 'Synthetic A', '4', '10', 'Off', 'meta'])
        code, titles = self.titles('run', root/'config.json', '--output', self.root/'segments', '--no-sample-label')
        self.assertEqual(code, 2)
        self.assertFalse(any('Synthetic A' in title for title in titles))
        fragment = json.loads((self.root/'segments/segment_001/preview_pair.json').read_text(encoding='utf-8'))
        self.assertIs(fragment['sample_label_in_title'], False)

    def test_iolite_export_help_names_the_provisional_run_directory(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as stopped:
            main(['iolite-export', '--help'])
        self.assertEqual(stopped.exception.code, 0)
        text = ' '.join(out.getvalue().split())
        self.assertIn('provisional_manifest.json with provisional_laser_log.csv', text)
        self.assertIn('correction_manifest.json with corrected_laser_log.csv', text)

    def test_figure_3_script_names_the_font_it_falls_back_to(self):
        # Arial is removed from the font list in memory, as on a system without it; no font file or cache is changed.
        package = Path(__file__).resolve().parents[1]
        script = package/'examples/R01/reproduce_figure3.py'
        code = ("import runpy, sys\nimport matplotlib\nmatplotlib.use('Agg')\nfrom matplotlib import font_manager\n"
                "manager = font_manager.fontManager\n"
                "manager.ttflist = [f for f in manager.ttflist if f.name.lower() != 'arial']\n"
                "manager._findfont_cached.cache_clear()\n"
                "sys.path.insert(0, sys.argv[1])\nsys.argv = [sys.argv[2], '--output', sys.argv[3]]\n"
                "runpy.run_path(sys.argv[0], run_name='__main__')\n")
        done = subprocess.run([sys.executable, '-B', '-X', 'utf8', '-c', code, str(script.parent), str(script), str(self.root/'figure')],
                              cwd=package, capture_output=True, text=True, encoding='utf-8')
        self.assertEqual(done.returncode, 0, done.stderr[-2000:])
        self.assertIn('Arial is not installed: the figure is drawn with DejaVu Sans.', done.stdout)
        self.assertIn('Audit matches expected_figure3.json', done.stdout)
        self.assertNotIn('findfont', done.stderr)
        self.assertTrue((self.root/'figure/Figure_3_R01_reproduced.png').is_file())


if __name__ == '__main__':
    unittest.main()
