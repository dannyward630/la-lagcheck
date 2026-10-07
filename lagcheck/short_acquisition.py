"""Conditional short-record analysis with the existing gates and fit unchanged."""
from .models import InputError


CHECKS = ('coverage', 'no_contradictory_high_quality_blocks', 'peak_not_at_search_boundary',
          'baseline_heldout', 'all_sensitivity_configs_estimable',
          'all_sensitivity_evidence_supported', 'sensitivity_envelope',
          'closure_converged', 'closure_coverage', 'closure_peak_not_at_boundary',
          'closure_no_contradictory_high_quality_blocks', 'closure_heldout', 'closure_residual')


def diagnose_short(ms, laser, settings, standard):
    """Use floor(N/12) blocks if needed; never reduce peak/evidence thresholds."""
    settings = dict(settings)
    switch = settings.pop('short_acquisition', laser.metadata.get('format') == 'new_wave_laser_log')
    if switch not in (True, False, 0, 1):
        raise InputError('short_acquisition must be true/false or 1/0')
    original, original_blocks = standard(ms, laser, settings)
    n = len(laser.lines)
    if not switch or n >= 180:
        return original, original_blocks
    requested = int(original['settings']['n_blocks'])
    available = n // 12
    used = min(requested, available)
    if used >= 4 and used != requested:
        report, blocks = standard(ms, laser, dict(settings, n_blocks=used))
    else:
        report, blocks = dict(original), original_blocks
    unresolved = None
    if report.get('reason') == 'initial_alignment_not_identifiable':
        candidate = report.get('alignment') or {}
        if candidate.get('offset_s') is not None and not candidate.get('boundary') and candidate.get('peak', 0) > .25 and used >= 4:
            unresolved = dict(candidate)
            report, blocks = standard(ms, laser, dict(settings, n_blocks=used, alignment_offset_s=candidate['offset_s']))
            report['alignment'] = dict(candidate, source='unverified_automatic_peak_candidate', used_for_conditional_inspection=True)
            report['checks']['initial_alignment_identifiable'] = False
    can_fit = report.get('baseline') is not None
    unavailable = [f'n_blocks={v}' for v in (8, 12, 15) if n < v * 12]
    execution = {}
    for key in CHECKS:
        if key not in report.get('checks', {}):
            execution[key] = dict(state='not_executed', reason=report.get('reason', 'no_fit'))
        elif key.startswith('all_sensitivity_') and unavailable:
            execution[key] = dict(state='incomplete', reason='too_few_intervals_for_' + ','.join(unavailable),
                                  result=report['checks'][key])
        else:
            execution[key] = dict(state='pass' if report['checks'][key] else 'fail', result=report['checks'][key])
    for variant in report.get('sensitivity', []):
        if variant['config'] in unavailable:
            variant['execution'] = 'not_executed'
            variant['unavailable_reason'] = 'fewer_than_12_intervals_per_block'
        elif variant.get('rate_ppm') is None:
            variant['execution'] = 'executed_no_estimate'
        else:
            variant['execution'] = 'pass' if variant['evidence_pass'] else 'fail'
    if unresolved is not None:
        execution['initial_alignment_identifiable'] = dict(state='fail', result=False, reason='coarse_peak_margin_below_0.01')
    report['standard_diagnostic'] = dict(status=original['status'], reason=original.get('reason'),
        requested_n_blocks=requested, rate_ppm=original.get('rate_ppm'), checks=original.get('checks', {}))
    report['short_acquisition'] = dict(enabled=True, policy='short_acquisition_v1',
        unverified_alignment_candidate=unresolved, interval_count=n, minimum_intervals_per_block=12, available_blocks=available,
        requested_blocks=requested, executed_blocks=len(blocks),
        usable_blocks=sum(bool(b.get('gate_pass') and b.get('fit_eligible')) for b in blocks),
        fitted_blocks=0 if not can_fit else report['baseline']['n_retained'],
        peak_min=original['settings']['peak_min'], margin_min=original['settings']['margin_min'],
        unavailable_sensitivity_configs=unavailable, check_execution=execution,
        explanation=('Short-record estimates retain the existing peak, margin, agreement, prediction and coverage gates. '
                     'The full 15-block sensitivity setting requires 180 intervals and was not executable. '
                     'Available checks and unexecuted checks are reported separately; a conditional estimate is not a verified correction.'))
    if report.get('rate_ppm') is not None:
        report['status'] = 'short_acquisition_conditional'
        failures = [k for k, v in execution.items() if v['state'] == 'fail']
        report['reason'] = 'Full sensitivity coverage unavailable: ' + ', '.join(unavailable)
        if failures:
            report['reason'] += '; executed checks failed: ' + ', '.join(failures)
    if unresolved is not None and report.get('rate_ppm') is None:
        report['reason'] += '; automatic alignment is unverified (coarse peak margin below 0.01)'
    return report, blocks
