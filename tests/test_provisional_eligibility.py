"""Eligibility of an estimate for provisional export.

A restriction on provisional output only. The rate, the checks, the status and
the exit code of a diagnosis are those of the estimator and are not changed.
"""
import contextlib
import copy
import csv
import io
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy.stats import t as student_t

from lagcheck.cli import run as _run, main, digest
from lagcheck.eligibility import assess, CONDITIONS, NOT_BLOCKED, STATUSES
from lagcheck.estimator import diagnose
from lagcheck.fallback import REQUIRED_CHECKS
from lagcheck.iolite import export_iolite
from lagcheck.laser import load_laser
from lagcheck.models import InputError, LaserData, SignalData
from lagcheck.ms import load_signal
from lagcheck.preview import render_preview_pair
from lagcheck.verification import verify_written_log
from examples.make_iolite_demo import make_demo
import test_cli

SPAN, TOL = 4000., .2
CLAIMS = ('stable', 'validated', 'verified', 'supported', 'accurate', 'reliable', 'passed')
INTERVAL_MS = {'path': 'signal.csv', 'format': 'csv', 'header_row': 0, 'encoding': 'utf-8', 'delimiter': ',',
               'time_kind': 'numeric', 'time_unit': 'ms', 'time_column': 'time_ms', 'signal_columns': ['element']}
INTERVAL_LASER = {'path': 'laser.csv', 'format': 'interval_csv', 'header_row': 0, 'encoding': 'utf-8', 'delimiter': ',',
                  'time_kind': 'numeric', 'time_unit': 's', 'start_column': 'start', 'end_column': 'end'}


def run(*args, **kwargs):
    kwargs.setdefault('preview', False)
    return _run(*args, **kwargs)


def diagnostic(rate=170., single_pass=None, error=1., blocks=10, alternatives=None, converged=True, history=None, **baseline):
    """A diagnostic with the fields that the eligibility reads, shaped as the estimator reports them."""
    single_pass = rate if single_pass is None else single_pass
    if alternatives is None:
        alternatives = {'n_blocks=8': single_pass+.5, 'n_blocks=12': single_pass-.5}
    variants = [dict(config=name, rate_ppm=value, evidence_pass=False) for name, value in alternatives.items()]
    variants.append(dict(config='closure', rate_ppm=rate if converged else None, evidence_pass=converged))
    if history is None:
        history = [dict(trial_ppm=single_pass, residual_ppm=rate-single_pass or .05), dict(trial_ppm=rate, residual_ppm=.05)]
    fit = dict(rate_ppm=single_pass, slope_se_ppm=error, n_retained=blocks, first_s=100., last_s=3900., max_residual_s=.01)
    return dict(rate_ppm=rate, baseline=dict(fit, **baseline), closure=dict(converged=converged, history=history), sensitivity=variants)


def wobbling(amplitude, rate=0.):
    """A record without drift whose block lags alternate by plus and minus ``amplitude`` seconds.

    The fitted rate then has no clear direction: 0.08 s gives a supported
    diagnosis, 0.12 s an iteration that does not converge, and 0.15 s also
    alternative settings whose rates differ in sign.
    """
    rng = np.random.default_rng(148)
    lines, cursor = [], 5.
    for _ in range(240):
        length = float(rng.uniform(5, 18))
        lines.append({'Start Time': cursor, 'End Time': cursor+length})
        cursor += length+float(rng.uniform(3, 10))
    t = np.arange(0, cursor+30, .1)
    signal = np.ones(len(t))*2.
    origin = lines[0]['Start Time']
    for i, row in enumerate(lines):
        lag = 3.+(amplitude if (i//24) % 2 == 0 else -amplitude)
        start = origin+(1+rate*1e-6)*(row['Start Time']-origin)+lag
        end = origin+(1+rate*1e-6)*(row['End Time']-origin)+lag
        signal += np.maximum(0, np.minimum(t+.05, end)-np.maximum(t-.05, start))/.1*1000
    signal += np.random.default_rng(149).normal(0, .05, len(t))
    return SignalData(t, signal), LaserData(lines)


def write_interval_case(root, ms, laser, **analysis):
    """Interval-CSV inputs and a configuration, as in test_cli (the readers rebase each time axis to its first value)."""
    root.mkdir(parents=True, exist_ok=True)
    with (root/'signal.csv').open('w', newline='') as f:
        w = csv.writer(f); w.writerow(['time_ms', 'element'])
        w.writerows(zip(ms.time*1000, ms.signal))
    with (root/'laser.csv').open('w', newline='') as f:
        w = csv.writer(f); w.writerow(['start', 'end'])
        w.writerows((row['Start Time'], row['End Time']) for row in laser.lines)
    config = {'ms': dict(INTERVAL_MS), 'laser': dict(INTERVAL_LASER),
              'analysis': dict({'tolerance_s': TOL, 'alignment_offset_s': 8., 'skip_record_edge_lines': False}, **analysis)}
    (root/'config.json').write_text(json.dumps(config), encoding='utf-8')
    return config


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


class RuleTests(unittest.TestCase):
    def assess(self, report, span=SPAN, tolerance=TOL, **kwargs):
        record = assess(report, span, tolerance, **kwargs)
        self.assertIn(record['status'], STATUSES)
        json.dumps(record, allow_nan=False)                 # the record is written to report.json as it is
        return record

    def test_clear_rate_is_not_blocked_and_the_numbers_of_each_condition_are_recorded(self):
        r = self.assess(diagnostic())
        self.assertEqual((r['status'], r['evidence_completeness']), ('not_blocked', 'complete'))
        self.assertEqual((r['triggered_conditions'], r['conditions_not_assessable'], r['missing_sensitivity_configs']), ([], [], []))
        self.assertEqual(r['candidate'], dict(source='acquisition_estimate', rate_ppm=170., rate_kind='iterated',
                                              single_pass_rate_ppm=170., single_pass_rate_se_ppm=1., single_pass_blocks=10))
        self.assertEqual((r['duration_s'], r['tolerance_s']), (SPAN, TOL))
        direction, sign, iteration = (r['conditions'][name] for name in CONDITIONS)
        quantile = float(student_t.ppf(.975, 8))
        self.assertEqual((direction['state'], direction['t_quantile'], direction['half_width_ppm']), ('not_triggered', quantile, quantile))
        self.assertEqual(direction['effect_s'], quantile*1e-6*SPAN)
        self.assertFalse(direction['includes_zero'])
        self.assertEqual((sign['state'], sign['min_rate_ppm'], sign['max_rate_ppm'], sign['signs_differ']), ('not_triggered', 169.5, 170.5, False))
        self.assertEqual(sign['effect_s'], 1.*1e-6*SPAN)
        self.assertEqual(iteration['state'], 'not_applicable')
        self.assertTrue(r['closure_checks_describe_candidate'])
        self.assertEqual(r['explanation'], NOT_BLOCKED)
        self.assertEqual(r['fit_coverage'], dict(coverage_fraction=None, first_block_mid_s=100., last_block_mid_s=3900.,
                                                 lines_left_out_at_start=None, lines_left_out_at_end=None))

    def test_conditional_short_record_is_assessed_and_lists_the_settings_it_cannot_execute(self):
        from test_new_wave_short import short_data
        ms, laser = short_data()                                     # 72 firing intervals: six blocks, three block counts not executable
        report, _ = diagnose(ms, laser, {'tolerance_s': TOL, 'alignment_offset_s': 3.})
        self.assertEqual(report['status'], 'short_acquisition_conditional')
        r = self.assess(report, span=laser.lines[-1]['End Time']-laser.lines[0]['Start Time'])
        self.assertEqual((r['status'], r['evidence_completeness'], r['missing_sensitivity_configs']),
                         ('not_blocked', 'incomplete', ['n_blocks=8', 'n_blocks=12', 'n_blocks=15']))
        self.assertEqual((r['candidate']['single_pass_blocks'], r['conditions']['rate_direction']['t_quantile']), (6, float(student_t.ppf(.975, 4))))
        self.assertEqual(r['fit_coverage']['coverage_fraction'], report['coverage_fraction'])

    def test_text_of_a_result_that_is_not_blocked_makes_no_claim_about_the_estimate(self):
        for report in (diagnostic(), diagnostic(alternatives={'n_blocks=8': 170.2, 'n_blocks=15': None}),
                       diagnostic(converged=False, history=[dict(trial_ppm=170., residual_ppm=.5), dict(trial_ppm=170.5, residual_ppm=.3)])):
            r = self.assess(report)
            self.assertEqual(r['status'], 'not_blocked')
            text = ' '.join([r['explanation'], *r['warnings'], *r['reasons']]).lower()
            for word in CLAIMS:
                self.assertNotIn(word, text)
        self.assertIn('minimum restriction on provisional output', NOT_BLOCKED)

    def test_negative_and_positive_rates_are_treated_alike(self):
        cases = [dict(rate=170., error=1.), dict(rate=4., error=30.),
                 dict(rate=40., single_pass=38., error=2., alternatives={'trim': -60., 'n_blocks=8': 41.})]
        for case in cases:
            mirrored = dict(case, rate=-case['rate'])
            for key in ('single_pass',):
                if key in case:
                    mirrored[key] = -case[key]
            if 'alternatives' in case:
                mirrored['alternatives'] = {k: -v for k, v in case['alternatives'].items()}
            a, b = self.assess(diagnostic(**case)), self.assess(diagnostic(**mirrored))
            self.assertEqual(a['status'], b['status'], case)
            for name in CONDITIONS:
                self.assertEqual(a['conditions'][name]['state'], b['conditions'][name]['state'], (case, name))
                self.assertEqual(a['conditions'][name].get('effect_s'), b['conditions'][name].get('effect_s'), (case, name))
        self.assertEqual(self.assess(diagnostic(rate=-170.))['status'], 'not_blocked')       # a stable negative rate is not an error
        self.assertEqual(self.assess(diagnostic(**cases[1]))['status'], 'withheld')
        self.assertEqual(self.assess(diagnostic(**cases[2]))['status'], 'withheld')

    def test_direction_needs_an_interval_that_includes_zero_and_an_effect_above_the_tolerance(self):
        near_zero = diagnostic(rate=.5, error=1.)                       # half-width 2.3 µs/s around a rate of 0.5 µs/s
        r = self.assess(near_zero, span=1000.)
        self.assertTrue(r['conditions']['rate_direction']['includes_zero'])
        self.assertEqual((r['conditions']['rate_direction']['state'], r['status']), ('not_triggered', 'not_blocked'))
        r = self.assess(near_zero, span=100000.)                        # the same uncertainty over a longer acquisition
        self.assertEqual((r['conditions']['rate_direction']['state'], r['status']), ('triggered', 'withheld'))
        self.assertEqual(r['triggered_conditions'], ['rate_direction'])
        self.assertIn('the direction of the single-pass rate is uncertain', r['reasons'][0])
        self.assertIn('above the tolerance of 0.2 s', r['reasons'][0])
        r = self.assess(diagnostic(rate=100., error=30.), span=100000.)     # a large effect, but the direction is clear
        self.assertFalse(r['conditions']['rate_direction']['includes_zero'])
        self.assertGreater(r['conditions']['rate_direction']['effect_s'], TOL)
        self.assertEqual(r['status'], 'not_blocked')

    def test_interval_that_just_touches_zero_includes_it(self):
        half_width = float(student_t.ppf(.975, 8))*2.
        touching = self.assess(diagnostic(rate=half_width, error=2.), span=1e6)
        self.assertTrue(touching['conditions']['rate_direction']['includes_zero'])
        self.assertEqual(touching['status'], 'withheld')
        beyond = self.assess(diagnostic(rate=math.nextafter(half_width, math.inf), error=2.), span=1e6)
        self.assertFalse(beyond['conditions']['rate_direction']['includes_zero'])
        self.assertEqual(beyond['conditions']['rate_direction']['state'], 'not_triggered')

    def test_effect_equal_to_the_tolerance_does_not_trigger_and_the_next_larger_effect_does(self):
        unconverged = dict(converged=False, alternatives={'n_blocks=8': 100.2, 'n_blocks=12': 100.1})
        cases = [('rate_direction', 'effect_s', diagnostic(rate=0., error=30.)),
                 ('sign_across_settings', 'effect_s', diagnostic(rate=5., alternatives={'trim': -20., 'n_blocks=8': 30.})),
                 ('unconverged_iteration', 'trial_difference_effect_s',
                  diagnostic(rate=100., history=[dict(trial_ppm=100., residual_ppm=60.), dict(trial_ppm=160., residual_ppm=.3)], **unconverged)),
                 ('unconverged_iteration', 'last_residual_effect_s',
                  diagnostic(rate=100., history=[dict(trial_ppm=100., residual_ppm=-70.)], **unconverged))]
        for name, key, report in cases:
            effect = self.assess(report)['conditions'][name][key]
            self.assertGreater(effect, 0, name)
            at_limit = self.assess(report, tolerance=effect)
            self.assertEqual((at_limit['conditions'][name]['state'], at_limit['status']), ('not_triggered', 'not_blocked'), (name, key))
            above = self.assess(report, tolerance=math.nextafter(effect, 0.))
            self.assertEqual((above['conditions'][name]['state'], above['status']), ('triggered', 'withheld'), (name, key))
            self.assertEqual(above['triggered_conditions'], [name])

    def test_zero_rate_with_zero_standard_error_is_an_exact_fit_and_not_a_failure(self):
        exact = diagnostic(rate=0., error=0., max_residual_s=0., alternatives={'n_blocks=8': 0., 'n_blocks=12': 0.})
        r = self.assess(exact)
        direction = r['conditions']['rate_direction']
        self.assertEqual((r['status'], direction['state'], direction['half_width_ppm'], direction['effect_s']), ('not_blocked', 'not_triggered', 0., 0.))
        self.assertTrue(direction['includes_zero'] and direction['zero_standard_error'])
        self.assertTrue(any('zero standard error does not include model error' in w for w in r['warnings']))
        # A zero that is not an exact fit over a time span cannot be used: the condition is not assessable.
        for change in (dict(max_residual_s=.01), dict(max_residual_s=None), dict(last_s=100.), dict(first_s=None)):
            r = self.assess(diagnostic(rate=0., error=0., **dict(dict(max_residual_s=0.), **change)))
            self.assertEqual((r['status'], r['conditions_not_assessable']), ('not_assessable', ['rate_direction']), change)
            self.assertIn('exact fit', r['conditions']['rate_direction']['not_available'][0])

    def test_missing_or_invalid_key_quantities_are_not_assessable_and_never_read_as_zero(self):
        invalid = [('rate_ppm', v) for v in (None, float('nan'), float('inf'), True, '170')]
        invalid += [('slope_se_ppm', v) for v in (None, float('nan'), float('inf'), -1., True)]
        invalid += [('n_retained', v) for v in (None, 2, 1, 0, True, 10.)]
        for key, value in invalid:
            r = assess(diagnostic(**{key: value}), SPAN, TOL)
            self.assertEqual(r['status'], 'not_assessable', (key, value))
            self.assertIn('rate_direction', r['conditions_not_assessable'], (key, value))
            self.assertIsNone(r['conditions']['rate_direction']['effect_s'])
            self.assertTrue(r['explanation'].startswith('A provisional log is not written from this estimate because its eligibility '
                                                        'could not be assessed'), r['explanation'])
        for span, tolerance in ((None, TOL), (0., TOL), (-1., TOL), (float('nan'), TOL), (float('inf'), TOL), (True, TOL),
                                (SPAN, None), (SPAN, 0.), (SPAN, float('nan')), (SPAN, True)):
            r = assess(diagnostic(), span, tolerance)
            self.assertEqual((r['status'], r['conditions_not_assessable']), ('not_assessable', ['rate_direction', 'sign_across_settings']),
                             (span, tolerance))
        r = assess(diagnostic(rate=float('inf')), SPAN, TOL)             # a rate that is not finite is not an estimate to write
        self.assertEqual(r['status'], 'not_assessable')
        self.assertIn('the estimated rate is not a finite number', r['reasons'])

    def test_condition_that_is_triggered_stands_when_another_cannot_be_assessed(self):
        r = self.assess(diagnostic(rate=40., slope_se_ppm=None, alternatives={'trim': -200., 'n_blocks=8': 41.}))
        self.assertEqual((r['status'], r['triggered_conditions'], r['conditions_not_assessable']),
                         ('withheld', ['sign_across_settings'], ['rate_direction']))
        self.assertEqual(len(r['reasons']), 2)

    def test_three_retained_blocks_leave_one_degree_of_freedom(self):
        r = self.assess(diagnostic(blocks=3))
        self.assertAlmostEqual(r['conditions']['rate_direction']['t_quantile'], 12.706204736, places=6)
        self.assertEqual(r['status'], 'not_blocked')

    def test_missing_alternative_is_listed_and_is_neither_zero_nor_a_reason_to_withhold(self):
        r = self.assess(diagnostic(alternatives={'n_blocks=8': 170.4, 'n_blocks=12': None, 'trim': None}))
        self.assertEqual((r['status'], r['evidence_completeness'], r['missing_sensitivity_configs']),
                         ('not_blocked', 'incomplete', ['n_blocks=12', 'trim']))
        sign = r['conditions']['sign_across_settings']
        self.assertEqual((sign['state'], sign['min_rate_ppm'], sign['signs_differ']), ('not_triggered', 170., False))
        warning = next(w for w in r['warnings'] if w.startswith('Evidence incomplete'))
        self.assertIn('n_blocks=12, trim', warning)
        self.assertIn('The sign of a missing estimate is unknown.', warning)

    def test_no_alternative_with_a_rate_is_not_assessable(self):
        none_estimable = diagnostic(alternatives={'n_blocks=8': None, 'trim': None}, converged=False,
                                    history=[dict(trial_ppm=170., residual_ppm=.5), dict(trial_ppm=170.5, residual_ppm=.3)])
        without_list = {k: v for k, v in diagnostic().items() if k != 'sensitivity'}
        for report in (none_estimable, without_list):
            r = self.assess(report)
            self.assertEqual((r['status'], r['evidence_completeness']), ('not_assessable', 'incomplete'))
            self.assertEqual(r['conditions_not_assessable'], ['sign_across_settings'])
            self.assertEqual(r['conditions']['sign_across_settings']['not_available'], ['a rate from at least one alternative setting'])

    def test_sign_condition_uses_every_finite_rate_and_needs_an_effect_above_the_tolerance(self):
        mixed = diagnostic(rate=40., alternatives={'trim': -60., 'n_blocks=8': 41.})
        r = self.assess(mixed, span=1000.)                    # rates of both signs, 0.1 s over the span
        sign = r['conditions']['sign_across_settings']
        self.assertEqual((sign['state'], sign['signs_differ'], sign['min_rate_ppm'], sign['max_rate_ppm']), ('not_triggered', True, -60., 41.))
        self.assertEqual(r['status'], 'not_blocked')
        r = self.assess(mixed)                                # the same rates, 0.4 s over the span
        self.assertEqual((r['status'], r['triggered_conditions']), ('withheld', ['sign_across_settings']))
        self.assertIn('differ in sign (-60 to 41 µs s⁻¹)', r['reasons'][0])
        # An alternative that failed its own evidence check still counts, and so does the rate that would be written.
        self.assertFalse(mixed['sensitivity'][0]['evidence_pass'])
        iterated_other_sign = self.assess(diagnostic(rate=-30., single_pass=30., alternatives={'n_blocks=8': 31., 'trim': 29.}))
        self.assertEqual(iterated_other_sign['conditions']['sign_across_settings']['min_rate_ppm'], -30.)
        self.assertEqual(iterated_other_sign['status'], 'withheld')
        all_negative = self.assess(diagnostic(rate=-300., alternatives={'trim': -100., 'n_blocks=8': -500.}))
        self.assertEqual((all_negative['status'], all_negative['conditions']['sign_across_settings']['signs_differ']), ('not_blocked', False))

    def test_standard_error_belongs_to_the_single_pass_rate_not_to_the_rate_after_iteration(self):
        # Single-pass rate near zero, rate after iteration far from it: the direction is judged on the single-pass fit.
        r = self.assess(diagnostic(rate=90., single_pass=4., error=30., alternatives={'n_blocks=8': 5., 'trim': 3.}))
        direction = r['conditions']['rate_direction']
        self.assertEqual((r['candidate']['rate_ppm'], r['candidate']['single_pass_rate_ppm'], r['candidate']['single_pass_rate_se_ppm']), (90., 4., 30.))
        self.assertEqual(direction['half_width_ppm'], float(student_t.ppf(.975, 8))*30.)
        self.assertTrue(direction['includes_zero'])
        self.assertEqual(r['triggered_conditions'], ['rate_direction'])
        # The other way round: a clear single-pass rate is not made uncertain by a small rate after iteration.
        r = self.assess(diagnostic(rate=4., single_pass=90., error=30., alternatives={'n_blocks=8': 91., 'trim': 89.}))
        self.assertFalse(r['conditions']['rate_direction']['includes_zero'])
        self.assertEqual(r['status'], 'not_blocked')

    def test_iteration_condition_applies_only_when_the_iteration_did_not_converge(self):
        wide_path = [dict(trial_ppm=100., residual_ppm=60.), dict(trial_ppm=160., residual_ppm=-.3), dict(trial_ppm=159.7, residual_ppm=.1)]
        converged = self.assess(diagnostic(rate=159.7, single_pass=100., alternatives={'n_blocks=8': 150.}, history=wide_path))
        self.assertEqual(converged['conditions']['unconverged_iteration'], dict(state='not_applicable', max_trial_difference_ppm=None,
                         trial_difference_effect_s=None, last_trial_ppm=None, last_residual_ppm=None, last_residual_effect_s=None, not_available=[]))
        self.assertEqual(converged['status'], 'not_blocked')         # a long path to convergence is not a reason to withhold
        self.assertTrue(converged['closure_checks_describe_candidate'])
        low = self.assess(diagnostic(rate=120., converged=False, alternatives={'n_blocks=8': 121., 'n_blocks=12': 119.},
                                     history=[dict(trial_ppm=120., residual_ppm=-2.), dict(trial_ppm=118., residual_ppm=.5)]), span=5000.)
        iteration = low['conditions']['unconverged_iteration']
        self.assertEqual((low['status'], low['evidence_completeness'], low['missing_sensitivity_configs']), ('not_blocked', 'incomplete', ['closure']))
        self.assertEqual(iteration['state'], 'not_triggered')
        self.assertAlmostEqual(iteration['trial_difference_effect_s'], 2.e-6*5000., places=12)
        self.assertAlmostEqual(iteration['last_residual_effect_s'], .5e-6*5000., places=12)
        self.assertEqual(iteration['last_trial_ppm'], 118.)
        self.assertFalse(low['closure_checks_describe_candidate'])
        self.assertTrue(any('are not evidence for the single-pass rate' in w and '118.000' in w for w in low['warnings']))
        high = self.assess(diagnostic(rate=150., converged=False, alternatives={'n_blocks=8': 151.},
                                      history=[dict(trial_ppm=150., residual_ppm=9.), dict(trial_ppm=159., residual_ppm=-2.)]), span=30000.)
        self.assertEqual((high['status'], high['triggered_conditions']), ('withheld', ['unconverged_iteration']))
        self.assertIn('differs from a trial rate of the iteration by up to 9 µs s⁻¹', high['reasons'][0])
        self.assertEqual(len(high['reasons']), 1)                    # the residual of 2 µs/s is 0.06 s and is recorded, not a reason
        self.assertLess(high['conditions']['unconverged_iteration']['last_residual_effect_s'], TOL)

    def test_unconverged_iteration_without_its_last_residual_is_not_assessable_unless_already_triggered(self):
        no_fit = [dict(trial_ppm=100., residual_ppm=1.), dict(trial_ppm=101., residual_ppm=None)]
        r = self.assess(diagnostic(rate=100., converged=False, alternatives={'n_blocks=8': 100.5}, history=no_fit))
        self.assertEqual((r['status'], r['conditions_not_assessable']), ('not_assessable', ['unconverged_iteration']))
        self.assertEqual(r['conditions']['unconverged_iteration']['not_available'], ['residual rate at the last trial'])
        far = [dict(trial_ppm=100., residual_ppm=90.), dict(trial_ppm=190., residual_ppm=None)]
        r = self.assess(diagnostic(rate=100., converged=False, alternatives={'n_blocks=8': 100.5}, history=far))
        self.assertEqual((r['status'], r['triggered_conditions']), ('withheld', ['unconverged_iteration']))
        unknown = diagnostic()
        unknown['closure'] = {}
        r = self.assess(unknown)
        self.assertEqual((r['status'], r['conditions_not_assessable']), ('not_assessable', ['unconverged_iteration']))

    def test_diagnosis_without_a_rate_has_no_estimate_to_assess(self):
        for report in ({'rate_ppm': None, 'reason': 'fewer_than_four_usable_blocks'}, {}, None):
            r = self.assess(report)
            self.assertEqual((r['status'], r['evidence_completeness'], r['triggered_conditions']), ('no_acquisition_rate', None, []))
            self.assertEqual({v['state'] for v in r['conditions'].values()}, {'not_applicable'})

    def test_relative_standard_error_above_one_half_alone_does_not_withhold(self):
        # A ratio of 0.6 with 0.06 s at stake is not withheld; a ratio of 0.45 with 0.4 s at stake is.
        small = self.assess(diagnostic(rate=10., error=6.))
        self.assertTrue(small['conditions']['rate_direction']['includes_zero'])
        self.assertEqual(small['status'], 'not_blocked')
        large = self.assess(diagnostic(rate=100., error=45., alternatives={'n_blocks=8': 101.}))
        self.assertEqual((large['status'], large['triggered_conditions']), ('withheld', ['rate_direction']))

    def test_standard_error_times_the_span_alone_does_not_decide(self):
        # A long record, standard error 6 µs/s: its effect over the span exceeds 0.2 s, but the direction is clear.
        long_record = self.assess(diagnostic(rate=151., single_pass=150., error=6., alternatives={'n_blocks=8': 145., 'trim': 153.}),
                                  span=40000.)
        self.assertGreater(6e-6*40000., TOL)
        self.assertEqual((long_record['status'], long_record['evidence_completeness']), ('not_blocked', 'complete'))
        # A short record, standard error 50 µs/s on seven blocks: its effect is below 0.2 s, but zero lies inside the half-width.
        short_record = self.assess(diagnostic(rate=90., error=50., blocks=7, converged=False, alternatives={'n_blocks=8': -700., 'trim': 100.},
                                              history=[dict(trial_ppm=90., residual_ppm=.4), dict(trial_ppm=90.4, residual_ppm=-.3)]), span=3000.)
        self.assertLess(50e-6*3000., TOL)
        self.assertEqual((short_record['status'], short_record['triggered_conditions']), ('withheld', ['rate_direction', 'sign_across_settings']))
        self.assertEqual(short_record['missing_sensitivity_configs'], ['closure'])


class Case(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def wobble(self, amplitude, name='case', **analysis):
        root = self.root/name
        write_interval_case(root, *wobbling(amplitude), **analysis)
        return root

    def interval(self, name='interval', step=0., declared=True, **analysis):
        """The interval-CSV fixture of test_cli: 240 lines and 100 µs/s, optionally with a 2 s step halfway."""
        root = self.root/name
        root.mkdir()
        config = test_cli.WorkflowTests.make_case(None, root, step)
        config['analysis'].update({'skip_record_edge_lines': False}, **analysis)
        if not declared:
            del config['analysis']['alignment_offset_s']
        (root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        return root

    def text(self, out):
        return (Path(out)/'report.md').read_text(encoding='utf-8').splitlines()

    def command(self, *args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = main([str(a) for a in args])
        return code, out.getvalue(), err.getvalue()

    def written_logs(self, out):
        return sorted(p.name for p in Path(out).iterdir() if p.name.endswith('laser_log.csv') or p.name.endswith('manifest.json'))

    def assert_no_claim(self, eligibility, text):
        """Neither the record nor its section of report.md describes an estimate that is not blocked as stable, validated or verified."""
        self.assertEqual(eligibility['status'], 'not_blocked')
        start = text.index('## Provisional export eligibility')
        end = next(i for i in range(start+1, len(text)) if text[i].startswith('## '))
        record = dict(eligibility, applied_mapping=dict(eligibility['applied_mapping'], file=None))
        for words in (' '.join(text[start:end]), json.dumps(record, ensure_ascii=False)):
            for word in CLAIMS:
                self.assertNotIn(word, words.lower(), word)


class WithheldEstimateTests(Case):
    def test_estimate_without_a_clear_direction_is_not_written_and_stays_in_the_report(self):
        root = self.wobble(.12)
        config = read(root/'config.json')
        expected, _ = diagnose(load_signal(root/'signal.csv', {k: v for k, v in config['ms'].items() if k != 'path'}),
                               load_laser(root/'laser.csv', {k: v for k, v in config['laser'].items() if k != 'path'}),
                               config['analysis'])                   # the diagnosis itself, without the command line
        code, r = run(root/'config.json', root/'out', True)
        self.assertEqual(code, 2)                                    # failed checks, not an input error
        for key in ('status', 'reason', 'rate_ppm', 'checks', 'failed_checks', 'baseline', 'sensitivity_envelope_s', 'coverage_fraction'):
            self.assertEqual(r.get(key), expected.get(key), key)
        self.assertEqual(r['status'], 'qc_failed')
        self.assertIsNotNone(r['rate_ppm'])
        self.assertEqual(self.written_logs(root/'out'), [])
        self.assertIsNone(r.get('provisional_file'))
        e = r['provisional_export_eligibility']
        self.assertEqual((e['status'], e['triggered_conditions'], e['evidence_completeness'], e['missing_sensitivity_configs']),
                         ('withheld', ['rate_direction'], 'incomplete', ['closure']))
        self.assertEqual((e['duration_s'], e['tolerance_s'], e['tolerance_request']), (r['duration_s'], TOL, dict(unit='s', value=TOL)))
        self.assertEqual(e['candidate']['single_pass_rate_se_ppm'], r['baseline']['slope_se_ppm'])
        self.assertEqual(e['fit_coverage'], dict(coverage_fraction=r['coverage_fraction'], first_block_mid_s=r['baseline']['first_s'],
                                                 last_block_mid_s=r['baseline']['last_s'], lines_left_out_at_start=0, lines_left_out_at_end=0))
        self.assertEqual(e['applied_mapping'], dict(kind='none', rate_ppm=None, rate_source=None, file=None, rejected_candidate=dict(
            rate_ppm=r['rate_ppm'], rate_source='acquisition_estimate', label='rejected candidate for inspection')))
        note = r['provisional_note']
        self.assertTrue(note.startswith('No provisional log was written from this acquisition estimate: '
                                        'the direction of the single-pass rate is uncertain'), note)
        self.assertTrue(note.endswith('The estimate and failed checks remain in the report. A documented reference rate may be supplied.'), note)
        self.assertNotIn('no rate could be estimated', note)         # a withheld estimate is not reported as a missing one
        self.assertEqual(e['message'], note)
        self.assertIn(note, r['notes'])
        text = self.text(root/'out')
        self.assertIn(f"Rate estimate (µs s⁻¹): {r['rate_ppm']}", text)
        self.assertIn('- closure_converged: not met', text)
        self.assertIn('Status: withheld; evidence incomplete', text)
        self.assertIn('- ' + note, text)
        self.assertIn(f"Rejected candidate for inspection: {r['rate_ppm']} µs s⁻¹ (acquisition_estimate). "
                      'It is not written to a corrected or provisional log and is not drawn as an applied correction.', text)
        self.assertEqual(read(root/'out/report.json')['provisional_export_eligibility'], e)
        code, printed, errors = self.command('run', root/'config.json', '--output', root/'cli', '--export', '--no-preview')
        summary = json.loads(printed)
        self.assertEqual((code, errors, summary['status'], summary['provisional_file'], summary['rate_ppm']), (2, '', 'qc_failed', None, r['rate_ppm']))
        self.assertIn(note, summary['notes'])

    def test_eligibility_is_recorded_without_export_and_with_provisional_output_disabled(self):
        root = self.wobble(.15)
        code, r = run(root/'config.json', root/'plain')
        e = r['provisional_export_eligibility']
        self.assertEqual((code, e['status'], e['triggered_conditions'], e['evidence_completeness']),
                         (2, 'withheld', ['rate_direction', 'sign_across_settings'], 'complete'))
        self.assertIsNone(e['message'])                               # no log was asked for, so none was withheld
        self.assertEqual(r['notes'], [])
        self.assertIn('differ in sign', e['explanation'])
        code, r = run(root/'config.json', root/'off', True, provisional=False)
        self.assertEqual((code, r['correction_note']), (2, 'No corrected log exported because timing checks did not pass'))
        self.assertIsNone(r['provisional_export_eligibility']['message'])
        self.assertEqual(self.written_logs(root/'off'), [])

    def test_reference_rate_replaces_a_withheld_estimate_and_takes_none_of_its_constants(self):
        root = self.wobble(.15)                                      # the iteration converges, so a fitted constant exists
        code, r = run(root/'config.json', root/'out', True, 'affine', '0', 'Known synthetic coordinate', reference_rate_ppm=170.9)
        self.assertEqual(code, 2)
        self.assertTrue(r['closure']['converged'] and r['closure']['residual_fit'])
        self.assertEqual((r['provisional_rate_source'], r['provisional_rate_ppm'], r['provisional_mode']), ('user_reference_rate', 170.9, 'rate_only'))
        plan, receipt = r['provisional_export_plan'], r['provisional_export_verification']
        self.assertEqual((plan['mode'], plan['offset_seconds'], plan['coordinate_reference']), ('rate_only', '0', None))
        self.assertEqual((receipt['rate_ppm'], receipt['offset_seconds'], receipt['alignment_offset_applied']), (170.9, '0', False))
        self.assertIn('a reference rate takes none; a rate-only provisional log was written instead', r['provisional_note'])
        self.assertIn(f"The rate estimated from this acquisition ({r['rate_ppm']:.3f} µs s⁻¹) was not written to a provisional log: ", r['provisional_note'])
        self.assertIn('This log applies the user-supplied reference rate (170.900 µs s⁻¹) as a rate_only mapping.', r['provisional_note'])
        use, e = r['reference_rate'], r['provisional_export_eligibility']
        self.assertEqual((use['used'], use['reason']), (True, 'acquisition_estimate_not_eligible'))
        self.assertNotIn('no rate could be estimated', use['note'])
        self.assertEqual(e['status'], 'withheld')
        self.assertTrue(e['message'].endswith('The provisional log applies the supplied reference rate instead.'))
        self.assertEqual((e['applied_mapping']['kind'], e['applied_mapping']['rate_source'], e['applied_mapping']['rejected_candidate']['rate_ppm']),
                         ('provisional_laser_log', 'user_reference_rate', r['rate_ppm']))
        manifest = read(root/'out/provisional_manifest.json')
        self.assertEqual((manifest['rate_source'], manifest['rate_ppm'], manifest['mode'], manifest['failure_stage']),
                         ('user_reference_rate', 170.9, 'rate_only', 'input_diagnostic'))
        self.assertEqual(manifest['provisional_export_eligibility']['status'], 'withheld')
        self.assertEqual(manifest['warnings'][0], e['message'])
        self.assertEqual(manifest['output_sha256'], digest(root/'out/provisional_laser_log.csv'))

    def test_reference_rate_in_the_three_situations(self):
        cases = [('none', self.interval('none', peak_min=.999), 'no_acquisition_rate', True, 'no_acquisition_estimate', 'user_reference_rate'),
                 ('withheld', self.wobble(.12, 'withheld'), 'withheld', True, 'acquisition_estimate_not_eligible', 'user_reference_rate'),
                 ('eligible', self.interval('eligible', step=2.), 'not_blocked', False, 'acquisition_estimate_available', 'acquisition_estimate')]
        for name, root, status, used, reason, source in cases:
            code, r = run(root/'config.json', root/'out', True, reference_rate_ppm=170.9)
            self.assertEqual(code, 2, name)
            self.assertEqual(r['provisional_export_eligibility']['status'], status, name)
            self.assertEqual((r['reference_rate']['used'], r['reference_rate']['reason']), (used, reason), name)
            self.assertEqual((r['provisional_rate_source'], r['provisional_export_eligibility']['applied_mapping']['rate_source']), (source, source), name)
            self.assertEqual(read(root/'out/provisional_manifest.json')['rate_source'], source, name)
        self.assertEqual(r['provisional_rate_ppm'], r['rate_ppm'])       # an eligible estimate is written, the reference rate is not

    def test_supported_diagnosis_keeps_its_verified_export_whatever_the_eligibility(self):
        root = self.wobble(.08)
        code, r = run(root/'config.json', root/'out', True)
        e = r['provisional_export_eligibility']
        self.assertEqual((code, r['status'], r['correction_exported'], e['status'], e['triggered_conditions']),
                         (0, 'supported', True, 'withheld', ['rate_direction']))
        self.assertEqual(self.written_logs(root/'out'), ['candidate_laser_log.csv', 'corrected_laser_log.csv', 'correction_manifest.json'])
        self.assertEqual((e['applied_mapping']['kind'], e['applied_mapping']['rate_ppm']), ('corrected_laser_log', r['rate_ppm']))
        self.assertIsNone(e['message'])
        text = self.text(root/'out')
        position = text.index('This result is supported. An export of it goes through the verified path, which this eligibility does not restrict.')
        # Being on the verified path is not being verified: the read-back of the written file decides that.
        self.assertEqual(text[position+1], 'A verified log still requires that the read-back of the written file passes.')
        # A supported result that is not blocked gets no such sentences, and its record makes no claim about the estimate.
        clear = self.interval('clear')
        code, r = run(clear/'config.json', clear/'out', True)
        self.assertEqual((code, r['status'], r['correction_exported']), (0, 'supported', True))
        text = self.text(clear/'out')
        self.assertFalse(any(line.startswith(('This result is supported.', 'A verified log still requires')) for line in text))
        self.assertIn(f"Mapping in use: corrected_laser_log.csv applies {r['rate_ppm']} µs s⁻¹ (acquisition_estimate).", text)
        self.assert_no_claim(r['provisional_export_eligibility'], text)

    def test_read_back_failure_of_a_candidate_that_is_not_eligible_publishes_nothing(self):
        root = self.wobble(.08)
        affine = (True, 'affine', '0', 'Known synthetic coordinate')
        code, r = run(root/'config.json', root/'out', *affine)
        self.assertEqual((code, r['input_diagnostic_status'], r['status']), (2, 'supported', 'export_qc_failed'))
        self.assertFalse(r['post_export_verification']['passed'])
        self.assertEqual(self.written_logs(root/'out'), ['candidate_laser_log.csv'])
        self.assertEqual(digest(root/'out/candidate_laser_log.csv'), r['export_verification']['output_sha256'])   # kept as checked
        self.assertTrue(r['provisional_note'].startswith('No provisional log was written from this acquisition estimate: '))
        self.assertTrue(r['correction_note'].startswith('No usable corrected file.'))
        self.assertIn(r['provisional_note'], r['notes'])
        self.assertEqual(r['provisional_export_eligibility']['applied_mapping']['kind'], 'none')
        # With a reference rate, that rate is written as its own rate-only log; the candidate is neither published nor corrected again.
        code, r = run(root/'config.json', root/'reference', *affine, reference_rate_ppm=170.9)
        self.assertEqual((code, r['status']), (2, 'export_qc_failed'))
        manifest = read(root/'reference/provisional_manifest.json')
        self.assertEqual((manifest['rate_source'], manifest['rate_ppm'], manifest['mode'], manifest['failure_stage']),
                         ('user_reference_rate', 170.9, 'rate_only', 'readback'))
        self.assertEqual(manifest['failed_checks'], r['post_export_verification']['failed_checks'])
        self.assertNotEqual(digest(root/'reference/provisional_laser_log.csv'), digest(root/'reference/candidate_laser_log.csv'))
        self.assertEqual(digest(root/'reference/candidate_laser_log.csv'), r['export_verification']['output_sha256'])
        self.assertIn('the input diagnostic passed, but the read-back of the candidate written from the acquisition estimate did not',
                      r['provisional_note'])
        self.assertEqual(r['provisional_readback']['status'] in ('passed', 'timing_checks_failed'), True)

    def test_eligible_candidate_is_still_published_unchanged_after_a_failed_read_back(self):
        import test_export_verification
        root = self.root/'affine'
        root.mkdir()
        test_export_verification.ExportVerificationTests.case(None, root)
        code, r = run(root/'config.json', root/'out', True, 'affine', '0', 'Known synthetic coordinate')
        self.assertEqual((code, r['status'], r['provisional_export_eligibility']['status']), (2, 'export_qc_failed', 'not_blocked'))
        self.assertEqual(digest(root/'out/provisional_laser_log.csv'), digest(root/'out/candidate_laser_log.csv'))    # no second correction
        manifest = read(root/'out/provisional_manifest.json')
        self.assertEqual((manifest['failure_stage'], manifest['rate_ppm']), ('readback', r['rate_ppm']))
        self.assertNotEqual(manifest['rate_ppm'], r['post_export_verification']['residual_rate_ppm'])    # the residual is not a new correction
        mapping = manifest['provisional_export_eligibility']['applied_mapping']
        self.assertEqual((mapping['kind'], mapping['rate_ppm'], mapping['rate_source'], mapping['file']),
                         ('provisional_laser_log', r['rate_ppm'], 'acquisition_estimate', r['provisional_file']))
        self.assertEqual(r['provisional_export_verification']['rate_ppm'], r['rate_ppm'])

    def test_changed_candidate_bytes_stay_an_execution_error_when_the_estimate_is_not_eligible(self):
        def change_after_readback(ms, laser, report, candidate, receipt, plan, analysis):
            result = verify_written_log(ms, laser, report, candidate, receipt, plan, analysis)
            self.assertFalse(result[0]['passed'])
            candidate.write_bytes(candidate.read_bytes() + b'\n')
            return result
        root = self.wobble(.08)
        with patch('lagcheck.cli.verify_written_log', side_effect=change_after_readback):
            with self.assertRaises(InputError):
                run(root/'config.json', root/'out', True, 'affine', '0', 'Known synthetic coordinate', reference_rate_ppm=170.9)
        self.assertTrue((root/'out/INCOMPLETE.json').exists())
        self.assertEqual(self.written_logs(root/'out'), ['candidate_laser_log.csv'])
        self.assertFalse((root/'out/report.json').exists())


class SelectedDiagnosticTests(Case):
    def diagnostic(self, reject=False, uncertain_retry=False):
        """diagnose() whose automatic result has no clear direction and fails the agreement checks.

        The retry at the input placement is the real one, unless ``uncertain_retry`` gives it the same uncertain fit.
        """
        def uncertain(report):
            report.update(rate_ppm=3.7)
            report['baseline'].update(rate_ppm=3.7, slope_se_ppm=30.)

        def implementation(ms, laser, analysis):
            report, blocks = diagnose(ms, laser, analysis)
            report = copy.deepcopy(report)
            if analysis.get('alignment_offset_s') is None:
                if report.get('baseline'):
                    uncertain(report)
                    report.update(status='qc_failed', reason='qc_failed', failed_checks=list(REQUIRED_CHECKS))
                    report['checks'] = dict(report['checks'], **{k: False for k in REQUIRED_CHECKS})
                self.strict = copy.deepcopy(report)
            else:
                if reject:
                    report['checks']['baseline_heldout'] = False
                if uncertain_retry:
                    uncertain(report)
                self.retry = copy.deepcopy(report)
            return report, blocks
        return implementation

    def test_eligibility_is_assessed_for_the_adopted_retry_not_for_the_strict_diagnostic(self):
        root = self.interval(declared=False)
        with patch('lagcheck.cli.diagnose', side_effect=self.diagnostic()):
            code, r = run(root/'config.json', root/'out', True)
        self.assertEqual(code, 2)
        for key in ('status', 'rate_ppm', 'checks', 'baseline'):
            self.assertEqual(r[key], self.strict[key], key)           # the strict fields at the top level are untouched
        self.assertTrue(r['alignment_fallback']['adopted'])
        e = r['provisional_export_eligibility']
        self.assertEqual(e['candidate'], dict(source='input_position_fallback', rate_ppm=self.retry['rate_ppm'], rate_kind='iterated',
                         single_pass_rate_ppm=self.retry['baseline']['rate_ppm'], single_pass_rate_se_ppm=self.retry['baseline']['slope_se_ppm'],
                         single_pass_blocks=self.retry['baseline']['n_retained']))
        self.assertEqual((e['status'], e['evidence_completeness']), ('not_blocked', 'complete'))
        self.assertEqual(assess(self.strict, e['duration_s'], TOL)['status'], 'withheld')    # the strict estimate alone would be withheld
        self.assertEqual((r['provisional_rate_ppm'], r['provisional_rate_source']), (self.retry['rate_ppm'], 'input_position_fallback'))
        self.assertEqual(r['distance_scale']['diagnostic'], 'input_position_fallback')
        self.assertEqual(read(root/'out/provisional_manifest.json')['provisional_export_eligibility']['candidate']['source'], 'input_position_fallback')

    def test_strict_diagnostic_is_assessed_when_the_retry_is_off_or_not_adopted(self):
        root = self.interval(declared=False)
        for name, kwargs, mock in (('off', dict(alignment_fallback=False), self.diagnostic()), ('rejected', {}, self.diagnostic(reject=True))):
            with patch('lagcheck.cli.diagnose', side_effect=mock):
                code, r = run(root/'config.json', root/name, True, **kwargs)
            e = r['provisional_export_eligibility']
            self.assertEqual((code, r['alignment_fallback']['adopted'], r['alignment_fallback']['attempted']), (2, False, name == 'rejected'), name)
            self.assertEqual((e['candidate']['source'], e['candidate']['single_pass_rate_ppm'], e['candidate']['single_pass_rate_se_ppm']),
                             ('acquisition_estimate', 3.7, 30.), name)
            self.assertEqual((e['status'], e['triggered_conditions']), ('withheld', ['rate_direction']), name)
            self.assertEqual(self.written_logs(root/name), [], name)
            self.assertIn('from this acquisition estimate', r['provisional_note'], name)
            self.assertEqual(r['rate_ppm'], 3.7, name)

    def test_adopted_retry_that_is_not_eligible_is_withheld_under_its_own_name(self):
        root = self.interval(declared=False)
        with patch('lagcheck.cli.diagnose', side_effect=self.diagnostic(uncertain_retry=True)):
            code, r = run(root/'config.json', root/'out', True)
        e = r['provisional_export_eligibility']
        self.assertEqual((code, r['alignment_fallback']['adopted'], e['candidate']['source'], e['status']), (2, True, 'input_position_fallback', 'withheld'))
        self.assertEqual(self.written_logs(root/'out'), [])
        self.assertTrue(r['provisional_note'].startswith('No provisional log was written from the estimate at the original input placement: '))
        self.assertEqual(e['applied_mapping']['rejected_candidate']['rate_source'], 'input_position_fallback')
        with patch('lagcheck.cli.diagnose', side_effect=self.diagnostic(uncertain_retry=True)):
            code, r = run(root/'config.json', root/'reference', True, reference_rate_ppm=170.9)
        self.assertEqual((r['provisional_rate_source'], r['provisional_failure_stage']), ('user_reference_rate', 'alignment_fallback'))
        self.assertIn('This log applies the user-supplied reference rate (170.900 µs s⁻¹)', r['provisional_note'])
        self.assertNotIn('This log applies the rate estimated at the original input placement', r['provisional_note'])
        self.assertEqual(r['provisional_export_plan']['offset_seconds'], '0')

    def test_fragments_are_assessed_separately_with_their_own_span_and_origin(self):
        root = make_demo(self.root/'demo')
        config = read(root/'config.json')
        del config['analysis']['alignment_offset_s']
        (root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        with (root/'laser.csv').open('a', encoding='utf-8', newline='') as f:
            w = csv.writer(f)
            w.writerows([['2026-01-01 02:00:00', 'Other', '2', '10', 'Off', 'metadata'],
                         ['2026-01-01 03:00:00', 'Synthetic A', '3', '10', 'On', 'later'], ['2026-01-01 03:00:01', '', '', '', 'Off', 'later']])
        with patch('lagcheck.cli.diagnose', side_effect=self.diagnostic()):
            code, r = run(root/'config.json', self.root/'segments', True, tolerance_s=.25)
        self.assertEqual(code, 2)
        first, later = (read(self.root/'segments'/name/'report.json') for name in ('segment_001', 'segment_002'))
        e = first['provisional_export_eligibility']
        self.assertEqual((e['candidate']['source'], e['status'], e['tolerance_s']), ('input_position_fallback', 'not_blocked', .25))
        self.assertEqual(e['duration_s'], first['alignment_fallback']['diagnostic']['duration_s'])
        self.assertEqual(first['tolerance_resolution']['request_source'], '--tolerance-s')
        self.assertEqual(first['tolerance_resolution']['overridden_configuration_value'], dict(unit='s', value=.2, source='analysis.tolerance_s'))
        self.assertTrue((self.root/'segments/segment_001/provisional_manifest.json').exists())
        self.assertEqual(later['alignment_fallback']['alignment_offset_s'], 10800.)       # the origin of the fragment is kept
        self.assertEqual(later['provisional_export_eligibility']['status'], 'no_acquisition_rate')
        self.assertFalse((self.root/'segments/segment_002/provisional_manifest.json').exists())   # the first fragment's result is not shared
        entries = read(self.root/'segments/segments.json')['segments']
        self.assertEqual([entry['provisional_export_eligibility'] for entry in entries], ['not_blocked', 'no_acquisition_rate'])
        # Each entry also gives the tolerance of its fragment in seconds. No speed converted a request given in seconds.
        self.assertEqual([(entry['resolved_tolerance_s'], entry['tolerance_conversion_speed_um_s']) for entry in entries], [(.25, None)]*2)


class CandidateEvidenceTests(Case):
    def test_unconverged_candidate_is_the_single_pass_rate_and_closure_checks_are_not_its_evidence(self):
        root = self.interval(step=2.)                                # the iteration ends at a trial rate that differs from the single-pass rate
        code, r = run(root/'config.json', root/'out', True)
        last = r['closure']['history'][-1]
        self.assertEqual((code, r['closure']['converged']), (2, False))
        self.assertGreater(abs(last['trial_ppm']-r['baseline']['rate_ppm']), 1.)
        self.assertTrue(r['checks']['closure_residual'] and r['checks']['closure_coverage'])    # checks met at the last trial rate
        e = r['provisional_export_eligibility']
        self.assertEqual((e['status'], e['evidence_completeness'], e['missing_sensitivity_configs']), ('not_blocked', 'incomplete', ['closure']))
        self.assertEqual((e['candidate']['rate_kind'], e['candidate']['rate_ppm']), ('single_pass', r['baseline']['rate_ppm']))
        self.assertIs(e['closure_checks_describe_candidate'], False)
        warning = next(w for w in e['warnings'] if w.startswith('The iteration did not converge.'))
        self.assertIn(f"evaluated at the last trial rate of the iteration ({last['trial_ppm']:.3f} µs s⁻¹) and are not evidence for "
                      'the single-pass rate.', warning)
        # The log and its fixed-mapping read-back use the single-pass rate, not the last trial rate.
        receipt, manifest = r['provisional_export_verification'], read(root/'out/provisional_manifest.json')
        self.assertEqual((r['provisional_rate_ppm'], receipt['rate_ppm'], manifest['rate_ppm']), (r['baseline']['rate_ppm'],)*3)
        self.assertEqual(manifest['output_sha256'], digest(root/'out/provisional_laser_log.csv'))
        self.assertEqual(read(root/'out/provisional_readback_report.json')['settings']['tolerance_s'], TOL)
        self.assertEqual(r['provisional_readback']['interpretation'], 'Re-diagnosis of the written provisional file with the mapping held fixed')
        mapping = e['applied_mapping']
        self.assertEqual((mapping['kind'], mapping['rate_ppm'], mapping['file']), ('provisional_laser_log', r['baseline']['rate_ppm'], r['provisional_file']))
        # The warnings travel with the log: report, manifest and report.md.
        self.assertEqual(manifest['warnings'], [NOT_BLOCKED, *e['warnings']])
        self.assertEqual(manifest['provisional_export_eligibility'], e)
        self.assertEqual(manifest['tolerance_resolution'], r['tolerance_resolution'])
        text = self.text(root/'out')
        self.assertIn('Status: not_blocked; evidence incomplete', text)
        for line in e['warnings']:
            self.assertIn(line, text)
        self.assertIn(NOT_BLOCKED, text)
        self.assertIn(f"Mapping in use: provisional_laser_log.csv applies {r['baseline']['rate_ppm']} µs s⁻¹ (acquisition_estimate).", text)
        self.assert_no_claim(e, text)
        # No comment line is added to the log: it has the header and one line for every record of the original.
        original = (root/'laser.csv').read_text(encoding='utf-8').splitlines()
        written = (root/'out/provisional_laser_log.csv').read_text(encoding='utf-8').splitlines()
        self.assertEqual((written[0], len(written)), (original[0], len(original)))

    def test_note_and_instruction_point_to_the_warnings_only_when_the_log_has_some(self):
        pointer = ' See also the warnings recorded in provisional_export_eligibility.warnings.'
        # A provisional log with warnings: the iteration did not converge, so one alternative setting gave no rate.
        root = self.interval(step=2.)
        code, r = run(root/'config.json', root/'out', True)
        manifest = read(root/'out/provisional_manifest.json')
        self.assertTrue(r['provisional_export_eligibility']['warnings'])
        self.assertEqual(manifest['provisional_export_eligibility']['warnings'], r['provisional_export_eligibility']['warnings'])
        for text in (r['provisional_note'], r['correction_note'], manifest['instruction']):
            self.assertTrue(text.endswith('do not apply the rate a second time.' + pointer), text)
            self.assertEqual(text.count('provisional_export_eligibility.warnings'), 1, text)
        self.assertEqual((r['correction_note'], manifest['instruction']), (r['provisional_note'],)*2)
        self.assertIn(r['correction_note'], self.text(root/'out'))
        # A provisional log without warnings keeps the three texts as they are without that sentence.
        plain = self.interval('plain', min_coverage=1.)
        code, r = run(plain/'config.json', plain/'out', True)
        manifest, e = read(plain/'out/provisional_manifest.json'), r['provisional_export_eligibility']
        self.assertEqual((code, e['status'], e['evidence_completeness'], e['warnings']), (2, 'not_blocked', 'complete', []))
        for text in (r['provisional_note'], r['correction_note'], manifest['instruction']):
            self.assertNotIn('provisional_export_eligibility.warnings', text)
            self.assertTrue(text.endswith('Determine the constant offset in the reconstruction importer; do not apply the rate a second time.'), text)
        # The other way to a provisional log, the checked candidate after a failed read-back, does the same.
        import test_export_verification
        for name, warnings in (('quiet', []), ('warned', ['A warning made for this test.'])):
            case = self.root/name
            case.mkdir()
            test_export_verification.ExportVerificationTests.case(None, case)
            with patch('lagcheck.cli.assess_eligibility', side_effect=lambda *args, **kwargs: dict(assess(*args, **kwargs), warnings=list(warnings))):
                code, r = run(case/'config.json', case/'out', True, 'affine', '0', 'Known synthetic coordinate')
            manifest = read(case/'out/provisional_manifest.json')
            self.assertEqual((code, r['provisional_failure_stage'], digest(case/'out/provisional_laser_log.csv')),
                             (2, 'readback', digest(case/'out/candidate_laser_log.csv')), name)
            for text in (r['provisional_note'], r['correction_note'], manifest['instruction']):
                self.assertTrue(text.startswith('Provisional, NOT verified: the input diagnostic passed, but the read-back of the written file did not'), text)
                self.assertEqual(text.endswith('state the failed checks with any result.' + pointer), bool(warnings), (name, text))
                self.assertEqual(text.endswith('state the failed checks with any result.'), not warnings, (name, text))

    def test_formatted_copy_of_a_provisional_log_keeps_its_status_and_warnings(self):
        inputs = make_demo(self.root/'inputs')
        config = read(inputs/'config.json')
        config['analysis']['n_blocks'] = 5                           # with five blocks the trimmed alternative has too few to give a rate
        (inputs/'config.json').write_text(json.dumps(config), encoding='utf-8')
        code, r = run(inputs/'config.json', self.root/'canonical', True)
        e = r['provisional_export_eligibility']
        self.assertEqual((code, e['status'], e['evidence_completeness'], e['missing_sensitivity_configs']), (2, 'not_blocked', 'incomplete', ['trim']))
        code, result = export_iolite(self.root/'canonical', self.root/'iolite')
        self.assertEqual((code, result['status']), (2, 'provisional_not_verified'))
        canonical, formatted = read(self.root/'canonical/provisional_manifest.json'), read(self.root/'iolite/provisional_manifest.json')
        self.assertEqual(canonical['warnings'], [NOT_BLOCKED, *e['warnings']])
        for key in ('provisional_export_eligibility', 'tolerance_resolution', 'warnings', 'rate_ppm', 'rate_source', 'failed_checks'):
            self.assertEqual(formatted[key], canonical[key], key)
        for key in ('provisional_export_eligibility', 'tolerance_resolution', 'warnings'):
            self.assertEqual(result[key], canonical[key], key)
        text = (self.root/'iolite/report.md').read_text(encoding='utf-8').splitlines()
        self.assertIn('Provisional export eligibility of the canonical run: not_blocked; evidence incomplete.', text)
        for line in canonical['warnings']:
            self.assertIn(line, text)


try:
    import matplotlib  # noqa: F401
    HAVE_MATPLOTLIB = True
except ImportError:
    HAVE_MATPLOTLIB = False


@unittest.skipUnless(HAVE_MATPLOTLIB, 'matplotlib is not installed')
class PreviewMappingTests(Case):
    def footer(self, *args, **kwargs):
        """Run and return (code, report, the footer text that was drawn)."""
        from matplotlib.figure import Figure
        with patch.object(Figure, 'supxlabel', autospec=True) as drawn:
            code, r = _run(*args, **kwargs)
        return code, r, [' '.join(call.args[1].split()) for call in drawn.call_args_list]

    def test_withheld_estimate_is_drawn_as_a_rejected_candidate_beside_the_uncorrected_view(self):
        root = self.wobble(.12)
        code, r, footers = self.footer(root/'config.json', root/'out', True)
        meta = read(root/'out/preview_pair.json')
        self.assertEqual((code, r['preview']['status']), (2, 'written'))
        self.assertEqual(meta['panels'], ['original / declared alignment', 'rejected candidate for inspection: acquisition_estimate, not applied'])
        self.assertEqual((meta['applied_rate_ppm'], meta['rate_source'], meta['mapping_provenance'], meta['first_line_lag_s']), (None, None, None, None))
        self.assertTrue(meta['rate_fitted'])                          # an estimate exists; it is not the applied mapping
        self.assertEqual(meta['rejected_candidate'], dict(rate_ppm=r['rate_ppm'], rate_source='acquisition_estimate', panel=meta['panels'][1]))
        starts = np.loadtxt(root/'laser.csv', delimiter=',', skiprows=1)[:, 0]
        self.assertTrue(np.allclose(meta['mapped_interval_starts_s'], starts-starts[0]+8., rtol=0, atol=1e-9))     # the declared alignment only
        self.assertEqual(meta['constant_lag_s'], 8.)
        self.assertFalse(any(w.startswith('No rate estimate') for w in meta['warnings']))
        self.assertIn('No rate correction is applied in this view. The rate estimated from this acquisition is drawn only as a '
                      'rejected candidate for inspection.', meta['warnings'])
        e = r['provisional_export_eligibility']
        self.assertEqual(meta['provisional_export_text'][:2], ['Provisional export eligibility: withheld, evidence incomplete. Tolerance 0.2 s.',
                                                               e['message']])
        self.assertEqual(meta['provisional_export_text'][2:], e['warnings'])
        for line in meta['provisional_export_text']:
            self.assertIn(' '.join(line.split()), footers[0])          # the text is on the figure, not only in the JSON
        self.assertEqual(meta['provisional_export_eligibility'], dict(status='withheld', evidence_completeness='incomplete',
                                                                      applied_mapping=e['applied_mapping']))
        # The view exists now, so the returned notes name it; report.json was written before the drawing.
        available = 'An uncorrected inspection view is available.'
        self.assertTrue(any(available in note for note in r['notes']))
        self.assertFalse(any(available in note for note in read(root/'out/report.json')['notes']))
        # A later preview of the same run directory reads the same record.
        again = render_preview_pair(root/'out', name='again')
        self.assertEqual((again['panels'], again['applied_rate_ppm']), (meta['panels'], None))

    def test_view_is_not_announced_when_the_drawing_fails(self):
        root = self.wobble(.12)
        with patch('lagcheck.preview.render_preview_pair', side_effect=ValueError('no drawing')), contextlib.redirect_stderr(io.StringIO()):
            code, r = _run(root/'config.json', root/'out', True)
        self.assertEqual(code, 1)
        self.assertNotEqual(r['preview']['status'], 'written')
        self.assertTrue(any(note.startswith('No provisional log was written') for note in r['notes']))
        self.assertFalse(any('An uncorrected inspection view is available.' in note for note in r['notes']))

    def test_figure_and_log_use_the_reference_rate_and_show_the_rejected_candidate_apart(self):
        root = self.wobble(.12)
        code, r, _ = self.footer(root/'config.json', root/'out', True, reference_rate_ppm=170.9)
        meta = read(root/'out/preview_pair.json')
        self.assertEqual(meta['panels'], ['one constant', 'correction: published_provisional_log',
                                          'rejected candidate for inspection: acquisition_estimate, not applied'])
        self.assertEqual((meta['applied_rate_ppm'], meta['mapping_provenance']['estimate_source'], meta['mapping_provenance']['sha256']),
                         (170.9, 'user_reference_rate', digest(r['provisional_file'])))
        self.assertEqual(meta['rejected_candidate']['rate_ppm'], r['rate_ppm'])
        # Without --export the same rate is displayed, as a rate only: the lag of the first line is the declared alignment.
        code, r, _ = self.footer(root/'config.json', root/'display', reference_rate_ppm=170.9)
        meta = read(root/'display/preview_pair.json')
        self.assertEqual((meta['panels'][1], meta['applied_rate_ppm'], meta['rate_source'], meta['first_line_lag_s']),
                         ('correction: user_reference_rate', 170.9, 'user_reference_rate', 8.))
        self.assertEqual(r['reference_rate']['reason'], 'export_not_requested')
        self.assertTrue(r['reference_rate']['note'].endswith('is not eligible for provisional export, so a preview of this run, '
                                                             'if one is drawn, displays the reference rate.'), r['reference_rate']['note'])
        # A run without a preview writes the same texts, so they do not say that a preview was drawn.
        code, quiet = run(root/'config.json', root/'quiet', reference_rate_ppm=170.9)
        self.assertFalse((root/'quiet/preview_pair.png').exists())
        self.assertEqual(quiet['reference_rate']['note'], r['reference_rate']['note'])
        self.assertIn('Mapping in use: no corrected or provisional log was written. A preview of this run, if one is drawn, displays '
                      '170.9 µs s⁻¹ (user_reference_rate).', self.text(root/'quiet'))
        run(root/'config.json', root/'none')
        self.assertIn('Mapping in use: no corrected or provisional log was written, and a preview of this run, if one is drawn, displays '
                      'no rate as applied.', self.text(root/'none'))

    def drawn(self, *args, **kwargs):
        """Run and return (code, report, the lines of the title, the footer) as they were put on the figure."""
        from matplotlib.figure import Figure
        with patch.object(Figure, 'suptitle', autospec=True) as titles, patch.object(Figure, 'supxlabel', autospec=True) as footers:
            code, r = _run(*args, **kwargs)
        return code, r, titles.call_args_list[0].args[1].splitlines(), ' '.join(footers.call_args_list[0].args[1].split())

    def test_figure_names_the_reference_rate_as_the_source_when_it_replaces_an_adopted_retry(self):
        root = self.interval(declared=False)
        with patch('lagcheck.cli.diagnose', side_effect=SelectedDiagnosticTests.diagnostic(self, uncertain_retry=True)):
            code, r, title, footer = self.drawn(root/'config.json', root/'out', True, reference_rate_ppm=170.9)
        meta, manifest = read(root/'out/preview_pair.json'), read(root/'out/provisional_manifest.json')
        written = digest(root/'out/provisional_laser_log.csv')
        self.assertEqual((code, r['alignment_fallback']['adopted'], r['provisional_export_eligibility']['status']), (2, True, 'withheld'))
        # The written file, its manifest, the record of the preview and the figure name one rate and one source.
        self.assertEqual((r['provisional_rate_ppm'], r['provisional_rate_source'], r['provisional_export_verification']['output_sha256']),
                         (170.9, 'user_reference_rate', written))
        self.assertEqual((manifest['rate_ppm'], manifest['rate_source'], manifest['output_sha256']), (170.9, 'user_reference_rate', written))
        self.assertEqual(meta['mapping_provenance'], dict(path=r['provisional_file'], sha256=written, mode='rate_only',
                                                          estimate_source='user_reference_rate'))
        self.assertEqual((meta['applied_rate_ppm'], meta['provisional_export_eligibility']['applied_mapping']['rate_source']),
                         (170.9, 'user_reference_rate'))
        self.assertEqual(meta['panels'], ['one constant', 'correction: published_provisional_log',
                                          'rejected candidate for inspection: input_position_fallback, not applied'])
        self.assertIn('The provisional log applies the supplied reference rate instead.', footer)
        # The read-back is that of the written file. This run wrote no candidate on the verified path.
        own = r['provisional_readback']
        self.assertEqual((meta['readback_status'], meta['readback_failed_checks'], meta['readback_record']),
                         (own['status'], own['failed_checks'], 'provisional_readback'))
        self.assertEqual((manifest['readback_status'], meta['candidate_readback_status'], meta['candidate_readback_failed_checks']),
                         (own['status'], None, []))
        self.assertEqual(title[1], f"Input: {r['input_diagnostic_status']}; read-back: {own['status']}; export: provisional_not_verified")

    def test_read_back_on_the_figure_is_that_of_the_displayed_log_and_that_of_the_candidate_is_stated_apart(self):
        root = self.wobble(.08)
        affine = (True, 'affine', '0', 'Known synthetic coordinate')
        _, first = run(root/'config.json', root/'first', *affine)
        # The estimate of this record is not eligible. Supplied as the reference rate, it is written as its own rate-only log.
        code, r, title, footer = self.drawn(root/'config.json', root/'out', *affine, reference_rate_ppm=first['rate_ppm'])
        meta, manifest = read(root/'out/preview_pair.json'), read(root/'out/provisional_manifest.json')
        post, own, written = r['post_export_verification'], r['provisional_readback'], digest(root/'out/provisional_laser_log.csv')
        self.assertEqual((code, r['status'], post['passed'], r['provisional_rate_source']), (2, 'export_qc_failed', False, 'user_reference_rate'))
        self.assertTrue(post['failed_checks'])
        self.assertNotEqual(written, digest(root/'out/candidate_laser_log.csv'))         # the displayed log is not the candidate
        # The read-back of the written file passed: the report, the manifest, the record of the preview and the title say so.
        self.assertEqual((own['status'], own['failed_checks'], manifest['readback_status']), ('passed', [], 'passed'))
        self.assertEqual((meta['readback_status'], meta['readback_failed_checks'], meta['readback_record']), ('passed', [], 'provisional_readback'))
        self.assertEqual((meta['mapping_provenance']['sha256'], meta['mapping_provenance']['estimate_source']), (written, 'user_reference_rate'))
        self.assertEqual((manifest['output_sha256'], manifest['rate_source']), (written, 'user_reference_rate'))
        self.assertEqual((meta['applied_rate_ppm'], manifest['rate_ppm']), (first['rate_ppm'],)*2)
        self.assertEqual(title[1:3], ['Input: supported; read-back of the displayed log: passed; export: provisional_not_verified',
                                     'Read-back of the candidate log (not the log displayed here): failed'])
        # The failed read-back of the candidate is kept, and it is named as that of the candidate.
        self.assertEqual((meta['candidate_readback_status'], meta['candidate_readback_failed_checks'], manifest['failed_checks']),
                         ('failed', post['failed_checks'], post['failed_checks']))
        lines = ['Failed read-back checks of the candidate log (written from the acquisition estimate; not the log displayed here): '
                 + ', '.join(post['failed_checks']),
                 'Read-back of the displayed provisional log (user_reference_rate), with its mapping held fixed: passed. This log stays provisional.']
        for line in lines:
            self.assertIn(line, meta['warnings'])
            self.assertIn(line, footer)
        self.assertFalse(any(w.startswith('Failed checks:') for w in meta['warnings']))
        # A read-back that passes does not make the log a verified one.
        self.assertEqual((meta['export_verification_status'], manifest['verified'], manifest['provisional'], r['correction_exported']),
                         ('provisional_not_verified', False, True, False))
        # report.md says the same: the failed stage and its checks are those of the candidate, the read-back is that of this file.
        text = (root/'out/report.md').read_text(encoding='utf-8')
        self.assertIn('\nFailed stage: readback of the candidate log written from the acquisition estimate (not this file)\n', text)
        self.assertIn('\nFailed read-back checks of that candidate log: ' + ', '.join(post['failed_checks']) + '\n', text)
        self.assertIn('\nFixed-mapping read-back of this provisional log: passed; residual rate of the written file', text)

    def test_log_with_the_bytes_of_the_candidate_and_another_mapping_states_its_own_read_back(self):
        from decimal import Decimal
        root = self.wobble(.08)
        _, first = run(root/'config.json', root/'first', True, 'affine', '0', 'Known synthetic coordinate')
        plan = first['export_plan']
        # With this output-clock origin the affine candidate is written with a constant of zero, so its bytes are those of a
        # rate-only log at the same rate. The two files still have different fixed mappings and different read-backs.
        origin = str(Decimal(plan['coordinate_reference']['target_origin_seconds'])-Decimal(plan['offset_seconds']))
        code, r, title, footer = self.drawn(root/'config.json', root/'out', True, 'affine', origin, 'Synthetic output-clock origin',
                                            reference_rate_ppm=first['rate_ppm'])
        meta, manifest = read(root/'out/preview_pair.json'), read(root/'out/provisional_manifest.json')
        post, own = r['post_export_verification'], r['provisional_readback']
        self.assertEqual(digest(root/'out/provisional_laser_log.csv'), digest(root/'out/candidate_laser_log.csv'))      # the same bytes
        self.assertEqual((code, r['export_plan']['mode'], r['provisional_export_plan']['mode'], r['provisional_rate_source']),
                         (2, 'affine', 'rate_only', 'user_reference_rate'))
        self.assertEqual((post['passed'], own['status'], own['failed_checks'], manifest['readback_status']), (False, 'passed', [], 'passed'))
        # The preview states the read-back of the log it displays, and that of the candidate apart.
        self.assertEqual((meta['readback_record'], meta['readback_status'], meta['readback_failed_checks']), ('provisional_readback', 'passed', []))
        self.assertEqual((meta['candidate_readback_status'], meta['candidate_readback_failed_checks']), ('failed', post['failed_checks']))
        self.assertEqual(meta['mapping_provenance']['estimate_source'], 'user_reference_rate')
        self.assertEqual(title[1:3], ['Input: supported; read-back of the displayed log: passed; export: provisional_not_verified',
                                     'Read-back of the candidate log (not the log displayed here): failed'])
        self.assertIn('Read-back of the displayed provisional log (user_reference_rate), with its mapping held fixed: passed. This log stays provisional.', footer)
        self.assertFalse(any(w.startswith('Failed checks:') for w in meta['warnings']))
        self.assertEqual((meta['export_verification_status'], manifest['verified'], r['correction_exported']), ('provisional_not_verified', False, False))

    def test_log_that_keeps_the_candidate_bytes_and_runs_without_such_a_log_state_the_read_back_as_before(self):
        import test_export_verification
        root = self.root/'affine'
        root.mkdir()
        test_export_verification.ExportVerificationTests.case(None, root)
        # An eligible candidate whose read-back failed is published unchanged, so its read-back is that of the displayed log.
        code, r, title, _ = self.drawn(root/'config.json', root/'out', True, 'affine', '0', 'Known synthetic coordinate')
        meta, post = read(root/'out/preview_pair.json'), r['post_export_verification']
        self.assertEqual(digest(root/'out/provisional_laser_log.csv'), digest(root/'out/candidate_laser_log.csv'))
        self.assertEqual((meta['readback_status'], meta['readback_failed_checks'], meta['readback_record']),
                         ('failed', post['failed_checks'], 'post_export_verification'))
        self.assertEqual((meta['candidate_readback_status'], meta['candidate_readback_failed_checks']), ('failed', post['failed_checks']))
        self.assertEqual(title[1], 'Input: supported; read-back: failed; export: provisional_not_verified')
        self.assertIn('Failed checks: ' + ', '.join(post['failed_checks']), meta['warnings'])
        self.assertEqual(meta['mapping_provenance']['estimate_source'], 'acquisition_estimate')
        text = (root/'out/report.md').read_text(encoding='utf-8')         # the wording of this case is as before
        self.assertIn('\nFailed stage: readback\nFailed checks: ' + ', '.join(post['failed_checks']) + '\nFixed-mapping read-back: ', text)
        # A verified log is the candidate, and its read-back passed.
        clear = self.interval('clear')
        code, r, title, _ = self.drawn(clear/'config.json', clear/'out', True)
        meta = read(clear/'out/preview_pair.json')
        self.assertEqual((meta['readback_status'], meta['readback_record'], meta['candidate_readback_status'], meta['mapping_provenance']['estimate_source']),
                         ('passed', 'post_export_verification', 'passed', 'acquisition_estimate'))
        self.assertEqual(title[1], 'Input: supported; read-back: passed; export: verified')
        # Without --export no read-back was asked for.
        code, r, title, _ = self.drawn(clear/'config.json', clear/'plain')
        meta = read(clear/'plain/preview_pair.json')
        self.assertEqual((meta['readback_status'], meta['readback_record'], meta['candidate_readback_status'], meta['candidate_readback_failed_checks']),
                         ('not_requested', None, None, []))
        self.assertEqual(title[1], 'Input: supported; read-back: not_requested; export: not_published')
        # A candidate that is not eligible is not published. The only read-back of that run is the one of the candidate.
        root = self.wobble(.08)
        code, r, title, _ = self.drawn(root/'config.json', root/'out', True, 'affine', '0', 'Known synthetic coordinate')
        meta = read(root/'out/preview_pair.json')
        self.assertEqual((meta['readback_status'], meta['readback_record'], meta['candidate_readback_status'], meta['export_verification_status']),
                         ('failed', 'post_export_verification', 'failed', 'not_published'))
        self.assertEqual(title[1], 'Input: supported; read-back: failed; export: not_published')
        # A reference-rate log of a record without a starting alignment: its read-back was not run.
        blind = self.interval('blind', declared=False, coarse_window_s=.01)
        code, r, title, _ = self.drawn(blind/'config.json', blind/'out', True, reference_rate_ppm=170.9, alignment_fallback=False)
        meta = read(blind/'out/preview_pair.json')
        self.assertEqual((r['provisional_readback']['status'], meta['readback_status'], meta['readback_record'], meta['candidate_readback_status']),
                         ('not_run', 'not_run', 'provisional_readback', None))
        self.assertEqual(read(blind/'out/provisional_manifest.json')['readback_status'], 'not_run')

    def test_report_without_the_mapping_record_is_drawn_as_before(self):
        from lagcheck.laser import load_laser
        from lagcheck.ms import load_signal
        from lagcheck.visualization import _map_lines
        root = self.wobble(.12)
        _, r = run(root/'config.json', root/'out')
        config = read(root/'out/config.json')
        laser = load_laser(config['laser']['path'], {k: v for k, v in config['laser'].items() if k != 'path'})
        ms = load_signal(config['ms']['path'], {k: v for k, v in config['ms'].items() if k != 'path'})
        recorded = _map_lines(config, r, laser, ms)
        self.assertEqual((recorded['rate'], recorded['source'], recorded['rejected']['rate']), (None, None, r['rate_ppm']))
        earlier = {k: v for k, v in r.items() if k != 'provisional_export_eligibility'}      # a report of an earlier version
        mapping = _map_lines(config, earlier, laser, ms)
        self.assertEqual((mapping['rate'], mapping['source'], mapping['rejected']), (r['rate_ppm'], 'acquisition_estimate', None))

    def test_supported_figure_has_no_added_text_and_other_results_state_tolerance_and_eligibility(self):
        root = self.interval()
        code, r, footers = self.footer(root/'config.json', root/'supported')
        meta = read(root/'supported/preview_pair.json')
        self.assertEqual((code, meta['provisional_export_text'], meta['rejected_candidate']), (0, [], None))
        self.assertNotIn('Provisional export eligibility', footers[0])
        root = self.interval('step', step=2.)
        code, r, footers = self.footer(root/'config.json', root/'out', True)
        meta = read(root/'out/preview_pair.json')
        self.assertEqual(meta['provisional_export_text'][0], 'Provisional export eligibility: not_blocked, evidence incomplete. Tolerance 0.2 s.')
        self.assertEqual(meta['provisional_export_text'][1], NOT_BLOCKED)
        self.assertIn('Evidence incomplete: no rate was obtained with the alternative setting(s) closure.', footers[0])
        self.assertEqual(meta['panels'], ['one constant', 'correction: published_provisional_log'])       # an eligible estimate is drawn as before


if __name__ == '__main__':
    unittest.main()
