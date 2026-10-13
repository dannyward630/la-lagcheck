"""Eligibility of an acquisition estimate for provisional export.

This is a minimum restriction on provisional output, separate from the timing
checks. It reads a completed diagnostic and changes none of its fields: the
rate, the checks and the status stay as the estimator reported them. The
status ``not_blocked`` means only that none of the three conditions below was
triggered.
"""
import math

STATUSES = ('not_blocked', 'withheld', 'not_assessable', 'no_acquisition_rate')
CONDITIONS = ('rate_direction', 'sign_across_settings', 'unconverged_iteration')
QUANTILE = .975
DEFINITIONS = ('T is the time from the first firing start to the last firing end of the selected acquisition. '
               'The screening half-width is the 0.975 quantile of Student\'s t with n - 2 degrees of freedom times the '
               'standard error of the single-pass slope, where n is the number of retained blocks. It is a scale of the '
               'conditional fit, not a confidence interval of the reported rate. A rate difference in µs s⁻¹ times T '
               'times 10⁻⁶ is the lag change in seconds that it amounts to over T. It is not a bound on the position error. '
               'T includes firing intervals at the ends of the record that are rewritten without having been in a fit, '
               'and the lag change over T is not an error bound for them. A condition is triggered only above the tolerance.')
NOT_BLOCKED = ('None of the three conditions that withhold a provisional log was triggered. This is a minimum restriction '
               'on provisional output. It changes no check, and the failed checks and missing estimates in this report '
               'still apply.')
LABELS = dict(rate_direction='direction of the single-pass rate',
              sign_across_settings='sign of the rate across the alternative settings',
              unconverged_iteration='change of the rate in an iteration that did not converge')


def _number(value):
    """A finite float, or None for anything else, booleans included."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def _positive(value):
    value = _number(value)
    return value if value is not None and value > 0 else None


def _over_span(rate_ppm, duration_s):
    """Lag change in seconds that a rate difference amounts to over the firing span, in the estimator's order of operations."""
    return rate_ppm*1e-6*duration_s


def _against(measured, limit):
    """Two numbers as text, with as few significant digits as still tell them apart (the rule of cli._against)."""
    for digits in range(4, 18):
        a, b = f'{measured:.{digits}g}', f'{limit:.{digits}g}'
        if a != b:
            break
    return a, b


def _above(effect, duration, tolerance):
    a, b = _against(effect, tolerance)
    return f'{a} s over the {duration:.0f} s from the first firing to the last, above the tolerance of {b} s'


def _exact_fit(baseline):
    """A zero standard error is accepted only for a fit with no residual over a nonzero time span."""
    first, last = _number(baseline.get('first_s')), _number(baseline.get('last_s'))
    return _number(baseline.get('max_residual_s')) == 0 and first is not None and last is not None and last > first


def _rate_direction(baseline, duration, tolerance):
    record = dict(state='not_assessable', t_quantile=None, half_width_ppm=None, includes_zero=None, effect_s=None,
                  zero_standard_error=None, not_available=[])
    rate, error, blocks = _number(baseline.get('rate_ppm')), _number(baseline.get('slope_se_ppm')), baseline.get('n_retained')
    if isinstance(blocks, bool) or not isinstance(blocks, int):
        blocks = None
    missing = [name for name, value in (('single-pass rate', rate), ('standard error of the single-pass slope', error),
                                        ('number of retained blocks', blocks), ('firing span', duration),
                                        ('tolerance', tolerance)) if value is None]
    if not missing and error < 0:
        missing.append('a standard error that is not negative')
    if not missing and blocks <= 2:
        missing.append('more than two retained blocks')
    if not missing and error == 0 and not _exact_fit(baseline):
        missing.append('evidence that the zero standard error is an exact fit')
    if missing:
        record['not_available'] = missing
        return record
    from scipy.stats import t as student_t
    quantile = float(student_t.ppf(QUANTILE, blocks-2))
    if not math.isfinite(quantile) or quantile <= 0:
        record['not_available'] = ['a finite t quantile']
        return record
    half_width = quantile*error
    effect = _over_span(half_width, duration)
    includes_zero = abs(rate) <= half_width      # an interval that just touches zero includes it
    record.update(state='triggered' if includes_zero and effect > tolerance else 'not_triggered', t_quantile=quantile,
                  half_width_ppm=half_width, includes_zero=includes_zero, effect_s=effect, zero_standard_error=error == 0)
    return record


def _sign_across_settings(rates, n_alternatives, duration, tolerance):
    record = dict(state='not_assessable', min_rate_ppm=None, max_rate_ppm=None, n_rates=len(rates), signs_differ=None,
                  effect_s=None, not_available=[])
    missing = [name for name, value in (('firing span', duration), ('tolerance', tolerance)) if value is None]
    if not n_alternatives:
        missing.append('a rate from at least one alternative setting')
    if missing:
        record['not_available'] = missing
        return record
    low, high = min(rates), max(rates)
    effect = _over_span(high-low, duration)
    differ = low < 0 < high
    record.update(state='triggered' if differ and effect > tolerance else 'not_triggered', min_rate_ppm=low, max_rate_ppm=high,
                  signs_differ=differ, effect_s=effect)
    return record


def _unconverged_iteration(closure, rate, duration, tolerance):
    record = dict(state='not_assessable', max_trial_difference_ppm=None, trial_difference_effect_s=None,
                  last_trial_ppm=None, last_residual_ppm=None, last_residual_effect_s=None, not_available=[])
    converged = closure.get('converged')
    if converged is True:
        record['state'] = 'not_applicable'
        return record
    if converged is not False:
        record['not_available'] = ['whether the iteration converged']
        return record
    history = closure.get('history') if isinstance(closure.get('history'), list) else []
    history = [step for step in history if isinstance(step, dict)]
    trials = [_number(step.get('trial_ppm')) for step in history]
    missing = [name for name, value in (('firing span', duration), ('tolerance', tolerance), ('rate to be written', rate))
               if value is None]
    if not missing:
        known = [trial for trial in trials if trial is not None]
        if known:
            difference = max(abs(rate-trial) for trial in known)
            record.update(max_trial_difference_ppm=difference, trial_difference_effect_s=_over_span(difference, duration))
        if not known or len(known) < len(trials):
            missing.append('every trial rate of the iteration')
        residual = _number(history[-1].get('residual_ppm')) if history else None
        if residual is None:
            missing.append('residual rate at the last trial')
        else:
            record.update(last_residual_ppm=residual, last_residual_effect_s=_over_span(abs(residual), duration))
        record['last_trial_ppm'] = trials[-1] if trials else None
    effects = [record['trial_difference_effect_s'], record['last_residual_effect_s']]
    if tolerance is not None and any(effect is not None and effect > tolerance for effect in effects):
        record['state'] = 'triggered'        # a known trigger stands even when the other quantity is not available
    elif not missing:
        record['state'] = 'not_triggered'
    record['not_available'] = missing
    return record


def _reasons(conditions, candidate, duration, tolerance):
    """One clause for each condition that was triggered, then one for each that could not be assessed."""
    reasons = []
    direction, sign, iteration = (conditions[name] for name in CONDITIONS)
    if direction['state'] == 'triggered':
        reasons.append(f"the direction of the single-pass rate is uncertain ({candidate['single_pass_rate_ppm']:.4g} µs s⁻¹ lies within "
                       f"the screening half-width of {direction['half_width_ppm']:.4g} µs s⁻¹ around zero, "
                       f"{candidate['single_pass_blocks']} blocks) and that half-width amounts to "
                       + _above(direction['effect_s'], duration, tolerance))
    if sign['state'] == 'triggered':
        reasons.append(f"the rates obtained with the alternative settings differ in sign ({sign['min_rate_ppm']:.4g} to "
                       f"{sign['max_rate_ppm']:.4g} µs s⁻¹) and that range amounts to " + _above(sign['effect_s'], duration, tolerance))
    if iteration['state'] == 'triggered':
        if iteration['trial_difference_effect_s'] is not None and iteration['trial_difference_effect_s'] > tolerance:
            reasons.append('the iteration did not converge and the rate to be written differs from a trial rate of the iteration by up to '
                           f"{iteration['max_trial_difference_ppm']:.4g} µs s⁻¹, which amounts to "
                           + _above(iteration['trial_difference_effect_s'], duration, tolerance))
        if iteration['last_residual_effect_s'] is not None and iteration['last_residual_effect_s'] > tolerance:
            reasons.append(f"the iteration did not converge and the residual rate at its last trial ({iteration['last_residual_ppm']:.4g} "
                           'µs s⁻¹) amounts to ' + _above(iteration['last_residual_effect_s'], duration, tolerance))
    for name in CONDITIONS:
        if conditions[name]['state'] == 'not_assessable':
            reasons.append(f"the {LABELS[name]} could not be assessed (not available: {', '.join(conditions[name]['not_available'])})")
    return reasons


def assess(fit_report, duration_s, tolerance_s, source='acquisition_estimate'):
    """Whether the rate of the selected diagnostic may be written to a provisional log.

    ``fit_report`` is the diagnostic an export would use: the adopted
    input-placement retry, otherwise the main diagnostic. ``duration_s`` is the
    time from the first firing start to the last firing end, and ``tolerance_s``
    the resolved tolerance of the run. The standard error is that of the
    single-pass fit and is used with the single-pass rate only. A quantity that
    is missing or not finite makes its condition ``not_assessable``; it is never
    read as zero and never as a condition that was not triggered.
    """
    report = fit_report if isinstance(fit_report, dict) else {}
    baseline = report.get('baseline') if isinstance(report.get('baseline'), dict) else {}
    closure = report.get('closure') if isinstance(report.get('closure'), dict) else {}
    duration, tolerance = _positive(duration_s), _positive(tolerance_s)
    rate = report.get('rate_ppm')
    candidate = dict(source=source, rate_ppm=rate,
                     rate_kind=None if rate is None else ('iterated' if closure.get('converged') else 'single_pass'),
                     single_pass_rate_ppm=baseline.get('rate_ppm'), single_pass_rate_se_ppm=baseline.get('slope_se_ppm'),
                     single_pass_blocks=baseline.get('n_retained'))
    # What the fit covers, shown beside T: T also spans firing intervals that no fit used.
    edge = report.get('record_edge_lines') if isinstance(report.get('record_edge_lines'), dict) else {}
    coverage = dict(coverage_fraction=report.get('coverage_fraction'), first_block_mid_s=baseline.get('first_s'),
                    last_block_mid_s=baseline.get('last_s'), lines_left_out_at_start=edge.get('skipped_at_start'),
                    lines_left_out_at_end=edge.get('skipped_at_end'))
    record = dict(status='no_acquisition_rate', evidence_completeness=None, candidate=candidate, duration_s=duration,
                  tolerance_s=tolerance, fit_coverage=coverage,
                  conditions={name: dict(state='not_applicable') for name in CONDITIONS},
                  triggered_conditions=[], conditions_not_assessable=[], missing_sensitivity_configs=[], reasons=[], warnings=[],
                  explanation='No rate was estimated from this acquisition, so no estimate was assessed.',
                  closure_checks_describe_candidate=None, definitions=DEFINITIONS)
    if rate is None:
        return record
    variants = report.get('sensitivity') if isinstance(report.get('sensitivity'), list) else []
    alternatives = [(str(v.get('config')), _number(v.get('rate_ppm'))) for v in variants if isinstance(v, dict)]
    missing = [name for name, value in alternatives if value is None]
    rates = [value for value in (_number(baseline.get('rate_ppm')), _number(rate)) if value is not None]
    rates += [value for _, value in alternatives if value is not None]
    conditions = dict(rate_direction=_rate_direction(baseline, duration, tolerance),
                      sign_across_settings=_sign_across_settings(rates, len(alternatives)-len(missing), duration, tolerance),
                      unconverged_iteration=_unconverged_iteration(closure, _number(rate), duration, tolerance))
    triggered = [name for name in CONDITIONS if conditions[name]['state'] == 'triggered']
    unassessed = [name for name in CONDITIONS if conditions[name]['state'] == 'not_assessable']
    reasons = _reasons(conditions, candidate, duration, tolerance)
    if _number(rate) is None:
        reasons.append('the estimated rate is not a finite number')
    incomplete = bool(missing or not alternatives)
    warnings = []
    if incomplete:
        warnings.append('Evidence incomplete: no rate was obtained with the alternative setting(s) '
                        + (', '.join(missing) or 'of the stability check') + '. The sign of a missing estimate is unknown. '
                        'The conditions were assessed on the available estimates only, and nothing is established about the missing ones.')
    unconverged = closure.get('converged') is False
    if unconverged:
        last = conditions['unconverged_iteration']['last_trial_ppm']
        warnings.append('The iteration did not converge. The candidate rate is the single-pass rate. The checks named closure_* '
                        'were evaluated at the last trial rate of the iteration'
                        + ('' if last is None else f' ({last:.3f} µs s⁻¹)') + ' and are not evidence for the single-pass rate.')
    if conditions['rate_direction'].get('zero_standard_error'):
        warnings.append('The standard error of the single-pass slope is zero: the retained blocks lie exactly on a line. '
                        'A zero standard error does not include model error.')
    if triggered:
        status, explanation = 'withheld', 'A provisional log is not written from this estimate: ' + '; '.join(reasons) + '.'
    elif unassessed or _number(rate) is None:
        status = 'not_assessable'
        explanation = ('A provisional log is not written from this estimate because its eligibility could not be assessed: '
                       + '; '.join(reasons) + '.')
    else:
        status, explanation = 'not_blocked', NOT_BLOCKED
    record.update(status=status, evidence_completeness='incomplete' if incomplete else 'complete',
                  conditions=conditions, triggered_conditions=triggered, conditions_not_assessable=unassessed,
                  missing_sensitivity_configs=missing, reasons=reasons, warnings=warnings, explanation=explanation,
                  closure_checks_describe_candidate=None if closure.get('converged') is None else not unconverged)
    return record
