"""A separate provisional estimate at the original input placement."""
import math


REQUIRED_CHECKS = ('no_contradictory_high_quality_blocks', 'baseline_heldout',
                   'closure_no_contradictory_high_quality_blocks', 'closure_heldout')
MANUAL_ADVICE = ('Manual next steps: declare a starting alignment only with an independent basis; '
                 'explicitly select a matching signal; or supply a reference rate from the same instrument pair. '
                 'See TROUBLESHOOTING.md.')


def try_alignment_fallback(ms, laser, strict, analysis, diagnose, *, enabled=True,
                           provisional=True, input_offset_s=0.):
    result = dict(attempted=False, adopted=False, status='not_attempted', triggered_reason=None,
                  alignment_offset_s=float(input_offset_s), required_checks=list(REQUIRED_CHECKS),
                  failed_adoption_checks=[], diagnostic=None, rate_source=None,
                  assumption='Original input placement; not an identified clock origin. '
                             'Correct pairing must lie within the configured block search range.')
    if not enabled or not provisional:
        result['note'] = 'Alignment fallback disabled by user.'
        return result, []
    if analysis.get('alignment_offset_s') is not None:
        result['note'] = 'User declared the starting alignment; no automatic fallback is attempted.'
        return result, []
    checks = strict.get('checks') or {}
    if strict.get('reason') == 'initial_alignment_not_identifiable' and strict.get('rate_ppm') is None:
        trigger = 'initial_alignment_not_identifiable'
    elif strict.get('rate_ppm') is not None and any(checks.get(k) is False for k in REQUIRED_CHECKS):
        trigger = 'automatic_estimate_failed_agreement_or_prediction'
    else:
        result['note'] = 'Strict diagnosis did not meet the alignment-fallback trigger.'
        return result, []
    candidate, blocks = diagnose(ms, laser, dict(analysis, alignment_offset_s=float(input_offset_s)))
    rate = candidate.get('rate_ppm')
    finite = isinstance(rate, (int, float)) and not isinstance(rate, bool) and math.isfinite(rate) and 1.+rate*1e-6 > 0
    failed = [k for k in REQUIRED_CHECKS if (candidate.get('checks') or {}).get(k) is not True]
    if not finite:
        failed.insert(0, 'finite_rate_with_positive_time_scale')
    adopted = not failed
    original_problem = ('Automatic start could not be identified. ' if trigger == 'initial_alignment_not_identifiable' else
                        'The automatic estimate failed block agreement or held-out prediction. ')
    note = original_problem + ('PROVISIONAL, NOT VERIFIED: a rate was estimated using the original input placement; '
            'this placement was not automatically identified. ' if adopted else
            'Input-placement fallback was attempted but not adopted: block agreement or held-out prediction '
            'did not pass, or a valid rate was unavailable. Details: ' + ', '.join(failed) + '. ')
    note += result['assumption'] + ' ' + MANUAL_ADVICE
    result.update(attempted=True, adopted=adopted,
                  status='provisional_adopted' if adopted else 'rejected', triggered_reason=trigger,
                  diagnostic=candidate, failed_adoption_checks=failed,
                  rate_source='input_position_fallback' if adopted else None, note=note)
    return result, blocks
