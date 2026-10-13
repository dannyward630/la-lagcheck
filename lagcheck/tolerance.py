"""Tolerance requests in seconds or micrometres, and distances equivalent to time differences.

This layer sits before and after the estimator. A request in micrometres is
resolved to seconds once, before the diagnosis, from the scan speed of the
valid firing intervals of the selected acquisition. The estimator, its checks
and its limits receive seconds only, and nothing here changes them.
"""
import math
import re

from .models import InputError
from .new_wave import FORMAT as NEW_WAVE_FORMAT
from .preview_geometry import SPEED_COLUMN, X_COLUMNS, Y_COLUMNS

# Settings of the analysis section that are used for reporting only and never reach the estimator.
REPORT_ONLY_KEYS = ('tolerance_um', 'spot_width_along_scan_um')
SPOT_SIZE_COLUMN, SPOT_TYPE_COLUMN, SPOT_ANGLE_COLUMN = 'Spot Size (um)', 'Spot Type', 'Spot Angle (deg)'
SECONDS_RULE = 'tolerance in seconds, used as given'
DISTANCE_RULE = ('tolerance in micrometres divided by the largest scan speed of the valid firing intervals; '
                 'resolved once, before the diagnosis')
NOT_CONVERTED = ('coverage fractions', 'correlation peaks and margins', 'numbers of settings that gave a rate',
                 'check results that are true or false', 'line-change intervals flagged by the laser-log screen')
_READ_ERRORS = (ValueError, TypeError, IndexError, KeyError, AttributeError)
# Where the speed in force at a firing start was written.
SPEED_ORIGINS = ('written_at_firing_start', 'inherited_from_record_above', 'not_recorded')
# Why the speed of a firing interval is not used, with the words of the message, in the order in which they are reported.
SPEED_NOT_USED = dict(
    not_recorded='no value at or above the firing start',
    not_a_positive_finite_number='magnitude not a positive finite number',
    changed_during_firing='changed during the firing',
    duration_not_positive='firing duration not positive',
    contradicted_by_positions='contradicted by the recorded positions and times by more than a factor of two',
    inherited_without_positions='taken from a record above the firing start, with no position of the interval to check it')
SCALE_NOTE = ('An equivalent distance is a time difference times the scan speed. It is a scale along a line scanned at that '
              'speed and does not establish the position error across a line change, a gap or a jump of the stage. '
              'A logged speed is a nominal value, not a measured one. A rate is multiplied by the firing span and by 10⁻⁶ '
              'before it is converted. A flagged line-change interval is a time deviation in the log, not a measured '
              'displacement. A spot width is a length scale only; it is not an image resolution.')


def estimator_settings(analysis):
    """The settings that the estimator accepts: keys used for reporting only are left out."""
    return {key: value for key, value in analysis.items() if key not in REPORT_ONLY_KEYS}


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise InputError(name + ' must be a positive finite number')
    return float(value)


def _number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def read_request(analysis, cli_s=None, cli_um=None, saved=None):
    """The tolerance request of this run: unit, value and where it was given.

    A source that gives both units is a conflict. One value on the command
    line replaces the tolerance of the configuration as a whole, unit included,
    and the replaced value is kept in the record. ``saved`` is the record that
    an earlier run wrote into its ``config.json``; its resolved seconds are
    reused and are not converted again.
    """
    analysis = analysis if isinstance(analysis, dict) else {}
    given = [key for key in ('tolerance_s', 'tolerance_um') if key in analysis]
    if len(given) == 2:
        raise InputError('analysis.tolerance_s and analysis.tolerance_um are both given; give the tolerance in one unit')
    if cli_s is not None and cli_um is not None:
        raise InputError('--tolerance-s and --tolerance-um are both given; give the tolerance in one unit')
    if saved is not None and not isinstance(saved, dict):
        raise InputError('tolerance_resolution must be the record written by an earlier run')
    configured = dict(unit=given[0][len('tolerance_'):], value=analysis[given[0]], source='analysis.' + given[0]) if given else None
    if cli_s is not None or cli_um is not None:
        unit, value = ('s', cli_s) if cli_s is not None else ('um', cli_um)
        return dict(unit=unit, value=_positive(value, '--tolerance-' + unit), source='--tolerance-' + unit,
                    overridden=configured, saved=None)
    request = dict(configured or dict(unit=None, value=None, source=None), overridden=None, saved=None)
    if request['unit'] == 'um':
        request['value'] = _positive(request['value'], 'analysis.tolerance_um')
    if saved is not None:
        if request['unit'] != 's' or saved.get('resolved_tolerance_s') != request['value']:
            raise InputError('tolerance_resolution records a tolerance that was resolved by an earlier run, and analysis.tolerance_s '
                             'does not repeat its resolved_tolerance_s; remove tolerance_resolution to give a new tolerance')
        request['saved'] = saved
    return request


def _inherited(rows, col):
    """Each cell of a column, a blank cell taking the last nonblank value above it (ESL convention)."""
    current, values = '', []
    for row in rows:
        cell = row[col].strip() if col < len(row) else ''
        current = cell or current
        values.append(current)
    return values


def _speed(text):
    """Magnitude of a recorded speed, or None when the cell is not a positive finite number."""
    try:
        value = abs(float(text))
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) and value > 0 else None


def _displacements(table, header):
    """Stage displacement (dx, dy) of each valid firing interval; None for an interval whose positions cannot be read.

    The result is None only when the log has no position columns at all. An interval
    with a blank, unreadable or non-finite coordinate has no position of its own, and
    what the positions of the other intervals show is kept.
    """
    xcol = next((c for c in X_COLUMNS if c in header), None)
    ycol = next((c for c in Y_COLUMNS if c in header), None)
    if xcol is None or ycol is None:
        return None
    ix, iy = header.index(xcol), header.index(ycol)
    moves = []
    for on, off in table['interval_rows']:
        try:
            sx, sy, ex, ey = (float(table['rows'][i][j]) for i in (on, off) for j in (ix, iy))
        except (ValueError, IndexError):
            moves.append(None)
            continue
        moves.append((ex-sx, ey-sy) if all(math.isfinite(v) for v in (sx, sy, ex, ey)) else None)
    return moves


def _log_speeds(table, header, lines, moves):
    """Speed in force at each firing start, where it was written and, when it is not used, the reason.

    A value written in the On record itself is explicit. A blank cell takes the
    last value above it. Such an inherited value can be a positioning speed, so it
    is used only when the positions and times of its own firing interval do not
    contradict it. An explicit value is used also when its interval has no position.
    """
    col = header.index(SPEED_COLUMN)
    rows, inherited = table['rows'], _inherited(table['rows'], col)
    speeds = []
    for i, (on, off) in enumerate(table['interval_rows']):
        own = rows[on][col].strip() if col < len(rows[on]) else ''
        text = own or inherited[on]
        origin = 'written_at_firing_start' if own else 'inherited_from_record_above' if text else 'not_recorded'
        value, reason = _speed(text), None
        duration = lines[i]['End Time']-lines[i]['Start Time']
        move = None if moves is None else moves[i]
        # A speed written in the Off record that ends the firing applies after it and is not compared.
        written = [rows[k][col].strip() for k in range(on, off) if col < len(rows[k]) and rows[k][col].strip()]
        if not text:
            reason = 'not_recorded'
        elif value is None:
            reason = 'not_a_positive_finite_number'
        elif any(_speed(cell) != value for cell in written):
            reason = 'changed_during_firing'
        elif not (math.isfinite(duration) and duration > 0):
            reason = 'duration_not_positive'
        elif move is not None:
            # A speed that the recorded positions and times of its interval contradict is not a firing speed,
            # for example a positioning speed that the firing start took from the record above it.
            if not .5*value*duration <= math.hypot(*move) <= 2.*value*duration:
                reason = 'contradicted_by_positions'
        elif origin != 'written_at_firing_start':
            reason = 'inherited_without_positions'
        speeds.append(dict(value=None if reason else value, origin=origin, not_used=reason))
    return speeds


def firing_speed(laser, declared=None):
    """Scan speed over the valid firing intervals of the selected acquisition, with its source.

    Only the speed in force at each firing start is read, so a speed written in
    an Off record that positions the stage is not taken as a firing speed.
    Several speeds are kept as a range. No speed is inferred from coordinates
    and durations. The positions of a firing interval are used for that
    interval only: to recognise a spot, and to set aside a recorded speed that
    they contradict by more than a factor of two. That factor is a coarse
    screen; a speed that passes it is not shown to be correct. A speed that the
    firing start took from a record above it is not used when its interval has
    no position to check it. A declared ``analysis.scan_speed_um_s`` is an
    explicit value of the user: it is used, and a difference from the log is
    recorded.
    """
    n = len(laser.lines)
    record = dict(status='not_computable', reason=None, source=None, n_firing_intervals=n, min_um_s=None, max_um_s=None,
                  values_um_s=[], declared_um_s=None, log_values_um_s=None, declared_differs_from_log=None,
                  declared_speed_accepted=True, log_speed_origins=None, log_speeds_not_used=None,
                  intervals_without_positions=None)
    fmt = (laser.metadata or {}).get('format')
    if fmt == NEW_WAVE_FORMAT:
        record.update(declared_speed_accepted=False,
                      reason='New Wave logs can carry a positioning speed at the first On record; no firing speed is taken from '
                             'this format and a declared speed is not applied to it')
        return record
    try:
        table = laser.table if isinstance(laser.table, dict) and fmt == 'laser_state_csv' else None
        header = [h.strip() for h in table['header']] if table else []
        if table and len(table['interval_rows']) != n:
            raise ValueError('the firing intervals and their records do not correspond')
        moves = _displacements(table, header) if table else None
        log = _log_speeds(table, header, laser.lines, moves) if table and SPEED_COLUMN in header else None
    except _READ_ERRORS:
        # Scan fields are read for reporting; a field that cannot be read never stops a run that gives seconds.
        # A declared speed is not applied to such a record, so the message of a distance request does not suggest one.
        record.update(declared_speed_accepted=False, reason='the scan fields of the laser log could not be read')
        return record
    if table:
        record['intervals_without_positions'] = n if moves is None else sum(move is None for move in moves)
    if log is not None:
        record.update(log_speed_origins={name: sum(speed['origin'] == name for speed in log) for name in SPEED_ORIGINS},
                      log_speeds_not_used={name: sum(speed['not_used'] == name for speed in log) for name in SPEED_NOT_USED})
    # One interval with known positions and no movement makes the record a spot record, whatever the other intervals hold.
    still = sum(move is not None and move[0] == 0 and move[1] == 0 for move in moves or [])
    if still:
        record.update(declared_speed_accepted=False,
                      reason=f'{still} of {n} firing intervals show no stage movement (spot record); a time difference has no '
                             'equivalent distance along a scan there')
        return record
    known = sorted({speed['value'] for speed in log or [] if speed['value'] is not None})
    if log is not None:
        record['log_values_um_s'] = known[:20]
    declared = _number(declared)
    if declared is not None and declared > 0:
        differs = any(not math.isclose(v, declared, rel_tol=1e-9) for v in known) if known else None
        record.update(status='credible', source='analysis.scan_speed_um_s (declared)', min_um_s=declared, max_um_s=declared,
                      values_um_s=[declared], declared_um_s=declared, declared_differs_from_log=differs)
        return record
    if log is None:
        record['reason'] = 'the laser log gives no scan speed and analysis.scan_speed_um_s is not set'
    elif any(speed['value'] is None for speed in log):
        counts = '; '.join(f'{count} {SPEED_NOT_USED[name]}' for name, count in record['log_speeds_not_used'].items() if count)
        record['reason'] = (f"{sum(speed['value'] is None for speed in log)} of {n} firing intervals have no usable scan speed "
                            f'in column {SPEED_COLUMN!r} ({counts}) and analysis.scan_speed_um_s is not set')
    else:
        record.update(status='credible', min_um_s=known[0], max_um_s=known[-1], values_um_s=record['log_values_um_s'],
                      source=f'laser log column {SPEED_COLUMN!r} at each firing start; a blank cell takes the last value above it')
    return record


def spot_width(laser, analysis):
    """Width of the spot along the scan direction with its source, or the reason it is not known.

    A declared ``analysis.spot_width_along_scan_um`` is used as given. From an
    ESL-style log a width is taken only when the shape, the size, the angle and
    the scan direction leave no choice: a circle, or a square with its sides
    along the stage axes that is scanned along an axis. The size must be a
    positive finite number and the angle of a square a finite number, and the
    three fields must stay the same from each firing start to its firing end.
    """
    record = dict(status='not_computable', width_um=None, source=None, reason=None)
    declared = analysis.get('spot_width_along_scan_um') if isinstance(analysis, dict) else None
    if declared is not None:
        record.update(status='known', width_um=_positive(declared, 'analysis.spot_width_along_scan_um'),
                      source='analysis.spot_width_along_scan_um (declared)')
        return record
    try:
        return _spot_width_from_log(laser, record)
    except _READ_ERRORS:
        record.update(status='not_computable', width_um=None, source=None, reason='the spot fields of the laser log could not be read')
        return record


def _shown(text, limit=40):
    """A cell as a message quotes it; a very long cell is cut."""
    return repr(text) if len(text) <= limit else f'{text[:limit]!r} and {len(text)-limit} more characters'


def _spot_width_from_log(laser, record):
    table = laser.table if isinstance(laser.table, dict) and (laser.metadata or {}).get('format') == 'laser_state_csv' else None
    header = [h.strip() for h in table['header']] if table else []
    if SPOT_SIZE_COLUMN not in header:
        record['reason'] = f'no column {SPOT_SIZE_COLUMN!r} is read from this laser log and analysis.spot_width_along_scan_um is not set'
        return record
    rows, intervals = table['rows'], table['interval_rows']
    columns = [header.index(name) if name in header else None for name in (SPOT_SIZE_COLUMN, SPOT_TYPE_COLUMN, SPOT_ANGLE_COLUMN)]
    inherited = [[''] * len(rows) if col is None else _inherited(rows, col) for col in columns]
    kinds = sorted(set(zip(*([values[on] for on, _ in intervals] for values in inherited))))
    if len(kinds) != 1:
        record['reason'] = 'the firing intervals record more than one spot size, type or angle'
        return record
    # The spot must stay as it is at the firing start until the firing ends. A cell written in the Off record
    # that ends the firing applies after it and is not compared.
    changed = sum(any(col is not None and col < len(rows[k]) and rows[k][col].strip() not in ('', values[on])
                      for col, values in zip(columns, inherited) for k in range(on+1, off))
                  for on, off in intervals)
    if changed:
        record['reason'] = (f'the spot size, type or angle changes during {changed} of {len(intervals)} firing intervals, '
                            'so the firing intervals have no single spot width')
        return record
    size, shape, angle = kinds[0]
    match = re.fullmatch(r'(\d+(?:\.\d+)?)(?:\s*[xX]\s*(\d+(?:\.\d+)?))?', size)
    sizes = [float(text) for text in match.groups() if text is not None] if match else []
    a, b = (sizes + [None])[:2] if sizes else (None, None)
    source = f'laser log columns {SPOT_SIZE_COLUMN!r}, {SPOT_TYPE_COLUMN!r} and {SPOT_ANGLE_COLUMN!r}'
    if not sizes or any(value <= 0 for value in sizes):
        record['reason'] = f'spot size {_shown(size)} is not read as one or two positive numbers'
    elif not all(math.isfinite(value) for value in sizes):
        # A long run of digits is read as a number and overflows; it is not a width.
        record['reason'] = f'spot size {_shown(size)} is not a finite number'
    elif not shape:
        record['reason'] = (f'the spot shape is not recorded ({SPOT_TYPE_COLUMN!r} is blank), so size {size!r} cannot be assigned '
                            'to the scan direction')
    elif shape.lower() == 'circle':
        if b is None or b == a:
            record.update(status='known', width_um=a, source=source + '; diameter of a circular spot')
        else:
            record['reason'] = f'a circular spot with two different sizes ({size!r}) is not interpreted'
    elif shape.lower() in ('square', 'rectangle'):
        moves = _displacements(table, header)
        try:
            turn = float(angle)
        except ValueError:
            turn = None
        if b is None or b != a:
            record['reason'] = f'which side of a spot of size {size!r} lies along the scan is not recorded'
        elif turn is None:
            record['reason'] = f'the spot angle is not recorded ({SPOT_ANGLE_COLUMN!r} is blank or not a number)'
        elif not math.isfinite(turn):
            # The remainder of a value that is not finite is not a number, and a comparison with it is false.
            record['reason'] = f'the spot angle {_shown(angle)} is not a finite number'
        elif min(turn % 90., 90.-turn % 90.) > 1e-9:
            record['reason'] = f'a square spot turned by {angle} degrees has no single width along the scan in this version'
        elif moves is None or any(move is None or min(abs(move[0]), abs(move[1])) > max(1e-6, 1e-6*max(abs(move[0]), abs(move[1])))
                                  for move in moves):
            record['reason'] = 'the scan direction is not recorded as a stage axis for every firing interval'
        else:
            record.update(status='known', width_um=a, source=source + '; side of a square spot aligned with the scan axis')
    else:
        record['reason'] = f'spot type {shape!r} is not interpreted'
    return record


def resolve(request, analysis, laser):
    """Estimator settings with the tolerance in seconds, the record of how it was obtained, and the speed evidence.

    Called once, before the diagnosis, on the original laser log. The
    read-back of a written log and a later export reuse
    ``resolved_tolerance_s``; they never derive a speed or a tolerance from
    the corrected log.
    """
    settings = estimator_settings(analysis)
    speed = firing_speed(laser, settings.get('scan_speed_um_s'))
    record = dict(requested_unit=request['unit'], requested_value=request['value'], request_source=request['source'],
                  overridden_configuration_value=request['overridden'], resolved_tolerance_s=None, rule=SECONDS_RULE,
                  speed_um_s=None, speed_source=None, reused_from_saved_configuration=False)
    if request['unit'] == 'um':
        if speed['status'] != 'credible':
            raise InputError('A tolerance in micrometres cannot be converted to seconds for this record: ' + speed['reason']
                             + '. Give the tolerance in seconds (analysis.tolerance_s or --tolerance-s)'
                             + (', or declare analysis.scan_speed_um_s' if speed['declared_speed_accepted'] else '') + '.')
        settings['tolerance_s'] = request['value']/speed['max_um_s']
        record.update(rule=DISTANCE_RULE, speed_um_s=speed['max_um_s'], speed_source=speed['source'])
    elif request['source'] == '--tolerance-s':
        settings['tolerance_s'] = request['value']
    if request.get('saved'):
        record = dict(request['saved'], reused_from_saved_configuration=True)
    else:
        record['resolved_tolerance_s'] = settings.get('tolerance_s')
    return settings, record, speed


def distance_scale(fit_report, blocks, tolerance_s, duration_s, speed, spot, eligibility=None, diagnostic='main'):
    """Second-valued quantities of the selected diagnostic with their equivalent distance along the scan.

    Only quantities that are time differences are converted. A rate enters as
    its lag change over the firing span. Coverage, correlation values, counts
    and true-or-false results have no equivalent distance and are left out.
    """
    fit = fit_report if isinstance(fit_report, dict) else {}
    low, high = (speed['min_um_s'], speed['max_um_s']) if speed['status'] == 'credible' else (None, None)
    width = _number(spot['width_um']) if spot['status'] == 'known' else None
    if width is not None and width <= 0:
        width = None        # a number of spot widths needs a positive finite width
    rows = []

    def add(quantity, seconds, basis='time difference', limit_s=None):
        seconds = _number(seconds)
        if seconds is None:
            return
        row = dict(quantity=quantity, seconds=seconds, limit_s=limit_s, basis=basis, um_min=None, um_max=None,
                   spot_widths_min=None, spot_widths_max=None)
        if low is not None:
            row.update(um_min=abs(seconds)*low, um_max=abs(seconds)*high)
            if width is not None:
                row.update(spot_widths_min=row['um_min']/width, spot_widths_max=row['um_max']/width)
        rows.append(row)

    tolerance, duration = _number(tolerance_s), _number(duration_s)
    baseline = fit.get('baseline') if isinstance(fit.get('baseline'), dict) else {}
    closure = fit.get('closure') if isinstance(fit.get('closure'), dict) else {}
    last_fit = closure.get('residual_fit') if isinstance(closure.get('residual_fit'), dict) else {}
    over_span = 'rate difference times the firing span'
    add('tolerance', tolerance)
    add('largest_heldout_error_single_pass', baseline.get('heldout_max_s'), limit_s=tolerance)
    add('largest_heldout_error_last_iteration', last_fit.get('heldout_max_s'), limit_s=tolerance)
    off_line = [abs(b['model_residual_s']) for b in blocks or [] if isinstance(b, dict)
                and b.get('block') in (fit.get('contradictory_high_quality_blocks') or []) and _number(b.get('model_residual_s')) is not None]
    add('largest_lag_difference_of_a_contradictory_block', max(off_line) if off_line else None, limit_s=tolerance)
    add('sensitivity_envelope', fit.get('sensitivity_envelope_s'), over_span, tolerance)
    residual = _number(closure.get('residual_ppm'))
    if residual is not None and duration is not None:
        add('closure_residual_over_firing_span', abs(residual)*1e-6*duration, over_span, tolerance)
    add('predicted_end_to_end_lag_change', fit.get('predicted_end_to_end_lag_change_s'), over_span)
    add('index_grid_lag_error_bound', fit.get('index_lag_error_bound_s'), limit_s=None if tolerance is None else tolerance/4)
    conditions = (eligibility or {}).get('conditions') or {}
    add('screening_half_width_over_firing_span', (conditions.get('rate_direction') or {}).get('effect_s'), over_span, tolerance)
    add('trial_difference_over_firing_span', (conditions.get('unconverged_iteration') or {}).get('trial_difference_effect_s'),
        over_span, tolerance)
    return dict(status='computed' if low is not None else 'not_computable', reason=speed['reason'], diagnostic=diagnostic,
                speed=speed, spot_width=spot, quantities=rows, not_converted=list(NOT_CONVERTED), note=SCALE_NOTE)
