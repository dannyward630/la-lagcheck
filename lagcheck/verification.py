"""Verify the actual exported bytes without silently realigning an affine file."""
from decimal import Decimal
import hashlib
import math
from pathlib import Path
import statistics

from .models import InputError, LaserData
from .laser import load_laser, _time
from .estimator import diagnose, SUPPORTED_STATUSES


def export_plan(ms, laser, report, options=None):
    options = dict(options or {})
    if set(options) - {'mode', 'ms_first_sample_at_laser_time', 'reference_note'}:
        raise InputError('Unknown export setting')
    mode = options.get('mode', 'rate_only')
    if mode not in ('rate_only', 'affine'):
        raise InputError('export.mode must be rate_only or affine')
    plan = dict(mode=mode, offset_seconds='0', coordinate_reference=None,
                initial_alignment_offset_s=report['alignment']['offset_s'])
    if mode == 'rate_only':
        if 'ms_first_sample_at_laser_time' in options:
            raise InputError('An explicit time reference requires export.mode=affine')
        return plan
    reference = options.get('ms_first_sample_at_laser_time')
    note = options.get('reference_note')
    if not isinstance(reference, str) or not reference.strip() or not isinstance(note, str) or not note.strip():
        raise InputError('Affine export needs ms_first_sample_at_laser_time and a reference_note explaining its source; no clock origin is guessed')
    target, template = _time(reference, laser.metadata['config'])
    original_template = laser.table['templates'][0][0]
    if template is not None and ((template.tzinfo is None) != (original_template.tzinfo is None)):
        raise InputError('The time reference and laser timestamps must have consistent timezone awareness')
    closed = report.get('closure', {})
    fit = closed.get('residual_fit')
    if not closed.get('converged') or not fit:
        raise InputError('Affine export requires a supported closure residual fit')
    retained = [b for b in closed['blocks'] if b['block'] in fit['retained']]
    anchor = statistics.median(float(b['t_mid']) for b in retained)
    residual_constant = fit['intercept_s'] + fit['rate_ppm'] * 1e-6 * (anchor-fit['intercept_time_s'])
    relative_offset = report['alignment']['offset_s'] + residual_constant
    # Convert from separately rebased coordinates to the declared output clock.
    offset = target - laser.table['origin'] + Decimal(str(relative_offset)) - Decimal(str(float(ms.time[0])))
    plan.update(offset_seconds=str(offset),
                coordinate_reference=dict(ms_first_sample_at_laser_time=reference,
                                          target_origin_seconds=str(target),
                                          reference_note=note,
                                          source='explicit_user_reference'),
                relative_alignment_offset_s=relative_offset,
                closure_constant_s=residual_constant,
                constant_anchor_time_s=anchor,
                constant_method='closure residual fitted lag at median retained block midpoint',
                import_instruction='Timestamps already include the fitted constant in the declared coordinate. Do not reapply the original offset or rate. Verify the receiving importer preserves this coordinate; zero additional offset is only a candidate, not a tested importer setting.')
    return plan


def verify_written_log(ms, laser, report, candidate, receipt, plan, analysis):
    candidate = Path(candidate)
    actual_hash = hashlib.sha256(candidate.read_bytes()).hexdigest()
    if actual_hash != receipt['output_sha256']:
        raise InputError('Candidate CSV changed after export; full verification cannot use these bytes')
    reloaded = load_laser(candidate, receipt['output_config'])
    settings = dict(analysis)
    if plan['mode'] == 'affine':
        target = Decimal(plan['coordinate_reference']['target_origin_seconds'])
        shift = float(reloaded.table['origin']-target) + float(ms.time[0])
        mapped = LaserData(lines=[dict(row, **{key: row[key]+shift for key in ('Start Time', 'End Time')})
                                 for row in reloaded.lines], metadata=reloaded.metadata, table=reloaded.table)
        # Fix zero after restoring the declared coordinate. A new automatic
        # alignment would hide an erroneous constant in the exported CSV.
        settings['alignment_offset_s'] = 0.0
    else:
        mapped = reloaded
        origin_change = float(reloaded.table['origin']-laser.table['origin'])
        settings['alignment_offset_s'] = report['alignment']['offset_s'] + origin_change
    rerun, blocks = diagnose(ms, mapped, settings)
    baseline_rate = (rerun.get('baseline') or {}).get('rate_ppm')
    reported_rate = rerun.get('rate_ppm')
    rates = (baseline_rate, reported_rate)
    finite = all(isinstance(x, (int, float)) and math.isfinite(x) for x in rates)
    rate_limit = report['settings']['closure_tolerance_ppm']
    tolerance = report['settings']['tolerance_s']
    duration = rerun.get('duration_s', math.inf)
    max_lag = max((abs(float(b['lag_s'])) for b in blocks if b.get('retained')), default=None)
    # A verified export needs an unqualified supported read-back. A read-back of a log written from a
    # qualified or failed diagnostic is informational; the qualified status then passes its own checks.
    accepted = ('supported',) if report.get('status') == 'supported' else SUPPORTED_STATUSES
    checks = dict(full_diagnostic_supported=rerun['status'] in accepted,
                  residual_rate_below_threshold=finite and max(abs(x) for x in rates)<rate_limit,
                  residual_drift_within_tolerance=finite and max(abs(x) for x in rates)*1e-6*duration<=tolerance,
                  coverage=rerun.get('coverage_fraction', 0)>=report['settings']['min_coverage'],
                  heldout=(rerun.get('baseline') or {}).get('heldout_max_s', math.inf)<=tolerance,
                  written_bytes_unchanged=hashlib.sha256(candidate.read_bytes()).hexdigest()==actual_hash)
    if plan['mode'] == 'affine':
        checks['fixed_reference_alignment_within_tolerance'] = max_lag is not None and max_lag<=tolerance
    summary = dict(passed=all(checks.values()), checks=checks,
                   failed_checks=[k for k,v in checks.items() if not v],
                   diagnostic_input_path=str(candidate.resolve()), diagnostic_input_sha256=actual_hash,
                   alignment_offset_s=settings['alignment_offset_s'], automatic_realignment_used=False,
                   residual_baseline_rate_ppm=baseline_rate, residual_rate_ppm=reported_rate,
                   residual_rate_limit_ppm=rate_limit, tolerance_s=tolerance,
                   coverage_fraction=rerun.get('coverage_fraction'),
                   heldout_max_s=(rerun.get('baseline') or {}).get('heldout_max_s'),
                   fixed_reference_max_lag_s=max_lag if plan['mode']=='affine' else None,
                   interpretation='Internal timing verification of this exported file, not independent spatial accuracy or a reconstruction-package round trip')
    return summary, rerun, blocks
