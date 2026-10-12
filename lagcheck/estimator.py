"""Block-lag estimation and acquisition-level timing diagnostics.

Supported means these specified internal checks passed, not independent accuracy.
This module never edits files and never substitutes a laboratory-wide rate.
"""
import contextlib
import io
import math
import numpy as np
from . import reference_core as core
from .models import InputError


def _native(value):
    if isinstance(value, dict):
        return {str(k): _native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_native(v) for v in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(value) else None
    return value


def _rescale(lines, rate, offset=0.0):
    origin = lines[0]['Start Time']
    factor = 1.0 + rate * 1e-6
    if not math.isfinite(factor) or factor <= 0:
        raise InputError('Estimated time scale must be finite and positive')
    return [dict(row, **{k: origin + factor * (row[k]-origin) + offset
                        for k in ('Start Time', 'End Time')}) for row in lines]


LINE_CHANGE_NEIGHBOURS = 5
LINE_CHANGE_SCALE_MULTIPLE = 10.0
LINE_CHANGE_DEFAULT_LIMIT_S = 1.0
BOUNDARY_DEFINITION = ('Search range cut short by the start or end of the signal record and peak within '
                       'two sampling intervals of that shortened boundary')
BOUNDARY_LIMITATION = ('The lag of a record-boundary block is measured near the boundary of its observable search range. '
                       'Agreement with the fitted line can show a lag inside the record, not a lag beyond the record boundary.')
SUPPORTED_STATUSES = ('supported', 'supported_excluding_record_boundary_blocks')


def screen_line_changes(lines, limit_s=1.0, neighbours=LINE_CHANGE_NEIGHBOURS,
                        scale_multiple=LINE_CHANGE_SCALE_MULTIPLE):
    """Compare each line-change interval of the laser log with its neighbours.

    A line-change interval is the logged start of firing on one interval minus
    the logged end of firing on the preceding one. An interval is flagged when
    it differs from the median of up to ``neighbours`` intervals on each side
    by more than ``limit_s`` and by more than ``scale_multiple`` times their
    scaled median absolute deviation. The second condition reduces flagging
    caused by ordinary variation in irregular sequences, such as spot
    acquisitions. A flagged interval may reflect a pause, a timestamp step or
    variation in the acquisition sequence; the log alone does not establish
    the cause. The block lags do not test an interval that lies inside an
    unusable block, so a flagged interval withholds supported status.
    """
    if isinstance(neighbours, bool) or not isinstance(neighbours, int) or neighbours < 1:
        raise InputError('neighbours must be a positive integer')
    for name, value in (('limit_s', limit_s), ('scale_multiple', scale_multiple)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise InputError(name + ' must be a positive finite number')
    start = np.array([float(row['Start Time']) for row in lines])
    end = np.array([float(row['End Time']) for row in lines])
    if not (np.isfinite(start).all() and np.isfinite(end).all()):
        raise InputError('Laser firing times must be finite')
    gaps = start[1:] - end[:-1]
    result = dict(method='neighbour_median', neighbours_each_side=int(neighbours), limit_s=float(limit_s),
                  default_limit_s=LINE_CHANGE_DEFAULT_LIMIT_S, scale_multiple=float(scale_multiple),
                  n_line_change_intervals=int(len(gaps)), evaluated=bool(len(gaps) >= 3), flagged=[],
                  n_flagged_at_default_limit=None,
                  median_interval_s=None, largest_unflagged_abs_deviation_s=None,
                  interpretation='A flag marks a line-change interval that is unusually long or short relative to its neighbours. '
                                 'It may reflect an acquisition pause, a timestamp step or variation in the acquisition sequence; '
                                 'the log alone does not establish the cause. Inspect the signal on either side and the acquisition '
                                 'design. Steps smaller than the limit are not detected by this screen.',
                  limitations='Intervals that recur by design, such as bracketing standards, are flagged. Groups of similarly '
                              'displaced adjacent intervals can escape detection because they shift the local reference. '
                              'Random scatter of several tenths of a second can produce a flag. '
                              'Intervals near the ends of a sequence have fewer neighbours.')
    if not result['evaluated']:
        return result
    medians = np.empty(len(gaps))
    scales = np.empty(len(gaps))
    for i in range(len(gaps)):
        near = np.r_[gaps[max(0, i-neighbours):i], gaps[i+1:i+1+neighbours]]
        medians[i] = np.median(near)
        scales[i] = 1.4826*np.median(np.abs(near-medians[i]))
    deviations = gaps - medians
    flagged = np.abs(deviations) > np.maximum(limit_s, scale_multiple*scales)
    # Recorded so that a raised limit stays visible in reports and manifests.
    result['n_flagged_at_default_limit'] = int((np.abs(deviations) > np.maximum(LINE_CHANGE_DEFAULT_LIMIT_S, scale_multiple*scales)).sum())
    result['flagged'] = [dict(after_interval=int(i+1), interval_s=float(gaps[i]),
                              neighbour_median_s=float(medians[i]), deviation_s=float(deviations[i]),
                              neighbour_scaled_mad_s=float(scales[i]), next_start_time_s=float(start[i+1]))
                         for i in np.flatnonzero(flagged)]
    result['median_interval_s'] = float(np.median(gaps))
    rest = np.abs(deviations[~flagged])
    result['largest_unflagged_abs_deviation_s'] = float(rest.max()) if len(rest) else None
    return result


def _fit(df, peak_min, margin_min, reject_tol=3.0, mode='worst'):
    if df is None or not len(df):
        return None
    good = df[(df.peak > peak_min) & (df.margin > margin_min)].copy()
    if 'fit_eligible' in good:
        good = good[good.fit_eligible]
    if mode == 'trim':
        good = good.iloc[1:-1]
    if len(good) < 4:
        return None
    x0 = float(good.t_mid.iloc[0])
    while True:
        x = good.t_mid.to_numpy()-x0
        y = good.lag_s.to_numpy()
        if np.ptp(x) <= 0:
            return None
        p = np.polyfit(x, y, 1)
        residual = y-np.polyval(p, x)
        if mode == 'all':
            keep = np.abs(residual) <= reject_tol
            if keep.sum() >= 4 and not keep.all():
                good = good.iloc[np.flatnonzero(keep)]
                mode = 'none'
                continue
            break
        if mode == 'none' or len(good) <= 4 or max(abs(residual)) <= reject_tol:
            break
        good = good.drop(good.index[int(np.argmax(abs(residual)))])
    errors = []
    for i in range(len(x)):
        mask = np.arange(len(x)) != i
        pred = np.polyval(np.polyfit(x[mask], y[mask], 1), x[i])
        errors.append(float(y[i]-pred))
    se = np.sqrt(np.sum(residual**2)/(len(x)-2)/np.sum((x-x.mean())**2))*1e6
    return dict(rate_ppm=float(p[0]*1e6), slope_se_ppm=float(se),
                intercept_s=float(p[1]), intercept_time_s=x0,
                retained=good.block.astype(int).tolist(), n_retained=len(good),
                first_s=float(good.t_mid.iloc[0]), last_s=float(good.t_mid.iloc[-1]),
                heldout_max_s=float(max(abs(np.array(errors)))),
                heldout_rms_s=float(np.sqrt(np.mean(np.array(errors)**2))),
                max_residual_s=float(max(abs(residual))))


# Largest lag-search half-width among the stability alternatives (max_shift = 30 s).
EDGE_GUARD_MIN_S = 30.0


def record_edge_lines(t, lines, guard, n_blocks, sides=('start', 'end')):
    """Count the lines at the named ends of the log that lie within ``guard`` seconds of the ends of the signal record.

    A block is matched to the signal by shifting its firing pattern, and the
    shift cannot run past the ends of the record. When the signal of the first
    line begins as the record begins, the first block cannot be shifted
    earlier, and its correlation peak is found at the edge of what can be
    searched. The diagnosis is then repeated without the lines that start
    within ``guard`` of the first recorded sample (or end within ``guard`` of
    the last), so that the first (or last) block has the signal of the
    neighbouring lines, and of the laser-off intervals between them, on both
    sides. ``guard`` is the largest search half-width used in the checks. The
    lines are chosen once, from the log at its starting alignment, and the
    same lines are left out in every fit of the repeated diagnosis. They are
    left out only when at least twelve lines per block remain.
    """
    guard = float(guard)
    lo, hi = float(t[0])+guard, float(t[-1])-guard
    head = 0
    while 'start' in sides and head < len(lines) and lines[head]['Start Time'] < lo:
        head += 1
    tail = 0
    while 'end' in sides and tail < len(lines)-head and lines[len(lines)-1-tail]['End Time'] > hi:
        tail += 1
    if (head or tail) and len(lines)-head-tail < 12*n_blocks:
        return dict(start=0, end=0, available_start=head, available_end=tail, applied=False, guard_s=guard)
    return dict(start=head, end=tail, available_start=head, available_end=tail, applied=bool(head or tail), guard_s=guard)


def _record_boundary_sides(df, closed_df, variants):
    """Ends of the record at which a block peak lies at a record-limited search boundary, in any fit."""
    sides = set()
    for frame in (df, closed_df):
        if frame is not None and 'record_boundary_block' in frame:
            sides.update(side for side in frame[frame.record_boundary_block].record_boundary_side if side)
    for variant in variants:
        sides.update(b['record_boundary_side'] for b in variant['blocks']
                     if b.get('record_boundary_block') and b.get('record_boundary_side'))
    return sides


def _measure(t, y, lines, settings, overrides=None):
    kw = {k: settings[k] for k in ('n_blocks', 'max_shift', 'peak_min', 'margin_min', 'reject_tol')}
    kw.update(overrides or {})
    with contextlib.redirect_stdout(io.StringIO()):
        _, df = core.analyse_loaded(t, y, lines, 'acquisition', **kw)
    if df is not None:
        dt=float((t[-1]-t[0])/(len(t)-1))
        pad=round(kw['max_shift']/dt)
        edges=np.linspace(0,len(lines),kw['n_blocks']+1).astype(int)
        lower=[]; upper=[]; observed=[]
        for row in df.to_dict(orient='records'):
            bi=int(row['block'])
            segment=lines[edges[bi]:edges[bi+1]]
            i0=int(np.searchsorted(t,segment[0]['Start Time']))
            i1=int(np.searchsorted(t,segment[-1]['End Time']))
            j0=max(0,i0-pad); j1=min(len(t),i1+pad)
            lower.append((j0-i0)*dt); upper.append((j1-i1)*dt)
            lag=row['lag_s']
            start=max(segment[0]['Start Time'],float(t[0]),float(t[0])-lag)
            end=min(segment[-1]['End Time'],float(t[-1]),float(t[-1])-lag)
            eligible=bool(np.isfinite(lag) and start<=segment[0]['Start Time']+1e-9
                          and end>=segment[-1]['End Time']-1e-9)
            observed.append(dict(observed_start_s=start,observed_end_s=max(start,end),
                                 observed_n_lines=sum(l['Start Time']>=start and l['End Time']<=end for l in segment),
                                 fit_eligible=eligible))
        df['search_min_s']=lower; df['search_max_s']=upper
        df['near_search_boundary']=(df.lag_s-df.search_min_s<=2*dt)|(df.search_max_s-df.lag_s<=2*dt)
        for key in ('observed_start_s','observed_end_s','observed_n_lines','fit_eligible'):
            df[key]=[row[key] for row in observed]
        # A record-boundary block has its search range cut short by the start
        # or end of the signal record and its peak within two sampling
        # intervals of that shortened boundary. When the record starts as the
        # signal of the first firing arrives, removing the fitted rate can
        # leave the first-block peak there; the search-edge criterion is then
        # not met.
        at_start=(df.search_min_s>-(pad-.5)*dt)&(df.lag_s-df.search_min_s<=2*dt)
        at_end=(df.search_max_s<(pad-.5)*dt)&(df.search_max_s-df.lag_s<=2*dt)
        boundary=np.isfinite(df.lag_s)&(at_start|at_end)
        df['record_boundary_block']=boundary
        df['record_boundary_side']=np.where(boundary&at_start,'start',np.where(boundary,'end',''))
        # Declared option: a fully observed record-boundary block leaves the
        # fit and the coverage, but its lag must still agree with the fitted
        # line (see _evidence). Blocks clipped by the record keep their rule.
        # Only a block that would otherwise enter the fit is excluded by the option.
        gate=(df.peak>kw['peak_min'])&(df.margin>kw['margin_min'])
        excluded=boundary&gate&df.fit_eligible.astype(bool)&bool(settings.get('exclude_record_boundary_blocks'))
        df['record_boundary_excluded']=excluded
        df.loc[excluded,'fit_eligible']=False
    return _fit(df, kw['peak_min'], kw['margin_min'], kw['reject_tol']), df


def _evidence(fit, df, lines, settings, require_coverage=True):
    """Use observed support, and retain contradictions in every fitted branch."""
    if fit is None:
        return dict(coverage_fraction=0.,coverage_required=require_coverage,peak_at_search_boundary=True,
                    contradictory_high_quality_blocks=[],evidence_pass=False)
    retained=df[df.block.isin(fit['retained'])]
    intervals=sorted(zip(retained.observed_start_s,retained.observed_end_s))
    total=0.; last=-np.inf
    for start,end in intervals:
        total+=max(0.,end-max(start,last));last=max(last,end)
    duration=lines[-1]['End Time']-lines[0]['Start Time']
    coverage=total/duration
    qualified=df[(df.peak>settings['peak_min']) & (df.margin>settings['margin_min'])]
    contradictory=[]
    for row in qualified.to_dict(orient='records'):
        if row['fit_eligible'] or row.get('record_boundary_excluded'):
            # An excluded record-boundary block is fully observed; it is held
            # to the same agreement tolerance at its midpoint.
            predicted=fit['intercept_s']+fit['rate_ppm']*1e-6*(row['t_mid']-fit['intercept_time_s'])
            discrepancy=abs(row['lag_s']-predicted)
        elif row['observed_end_s']>row['observed_start_s']:
            # A clipped block cannot define a reliable nominal midpoint, but
            # its observed high-quality lag can still contradict the model.
            # Allow every lag the fitted line reaches over the observed span.
            predictions=[fit['intercept_s']+fit['rate_ppm']*1e-6*(x-fit['intercept_time_s'])
                         for x in (row['observed_start_s'],row['observed_end_s'])]
            discrepancy=max(0.,min(predictions)-row['lag_s'],row['lag_s']-max(predictions))
        else:
            continue
        if discrepancy>settings['tolerance_s']:
            contradictory.append(int(row['block']))
    boundary=bool(retained.near_search_boundary.any())
    return dict(coverage_fraction=coverage,coverage_required=require_coverage,peak_at_search_boundary=boundary,
                contradictory_high_quality_blocks=contradictory,
                evidence_pass=bool((not require_coverage or coverage>=settings['min_coverage']) and not boundary and not contradictory
                                   and fit['heldout_max_s']<=settings['tolerance_s']))


def _boundary_rows(frame, model):
    if frame is None or 'record_boundary_block' not in frame:
        return []
    retained = None if model is None else frame[frame.block.isin(model['retained'])]
    median_peak = None if retained is None or not len(retained) else float(retained.peak.median())
    rows = []
    for row in frame[frame.record_boundary_block].to_dict(orient='records'):
        start = row['record_boundary_side'] == 'start'
        distance = row['lag_s']-row['search_min_s'] if start else row['search_max_s']-row['lag_s']
        predicted = None if model is None else model['intercept_s']+model['rate_ppm']*1e-6*(row['t_mid']-model['intercept_time_s'])
        rows.append(dict(block=int(row['block']), side=row['record_boundary_side'], lag_s=row['lag_s'], peak=row['peak'],
                         median_peak_of_retained_blocks=median_peak, distance_to_record_boundary_s=distance,
                         lag_minus_fitted_line_s=None if predicted is None else row['lag_s']-predicted,
                         excluded_from_fit_by_option=bool(row['record_boundary_excluded'])))
    return rows


def _boundary_report(settings, df, fit, closed_df, closed, variants):
    """List record-boundary blocks of the reported fits and of the stability alternatives, under either policy."""
    report = dict(policy='exclude' if settings['exclude_record_boundary_blocks'] else 'check',
                  baseline=_boundary_rows(df, fit), closure=_boundary_rows(closed_df, closed),
                  definition=BOUNDARY_DEFINITION, limitation=BOUNDARY_LIMITATION)
    report['excluded_blocks'] = sorted({r['block'] for r in report['baseline']+report['closure']
                                        if r['excluded_from_fit_by_option']})
    listed, excluded = {}, {}
    for variant in variants:
        if variant['config'] in ('all', 'trim', 'closure'):
            continue        # these reuse the blocks already listed under baseline or closure
        here = [int(b['block']) for b in variant['blocks'] if b.get('record_boundary_block')]
        gone = [int(b['block']) for b in variant['blocks'] if b.get('record_boundary_excluded')]
        if here:
            listed[variant['config']] = here
        if gone:
            excluded[variant['config']] = gone
    report['in_sensitivity'] = listed
    report['excluded_in_sensitivity'] = excluded
    return report


def _coarse(t, y, lines, window):
    log_on, wave = core.on_sequences(t, y, lines)
    a, b = log_on-log_on.mean(), wave-wave.mean()
    den = np.linalg.norm(a)*np.linalg.norm(b)
    if not den > 0:
        return dict(ok=False, reason='no_temporal_contrast')
    size = 1 << (2*len(t)-1).bit_length()
    corr = np.fft.irfft(np.conj(np.fft.rfft(a, size))*np.fft.rfft(b, size), size)/den
    dt = float((t[-1]-t[0])/(len(t)-1))
    count = min(int(window/dt), len(t)-2)
    if count < 2:
        return dict(ok=False, reason='coarse_window_too_small')
    lag = np.arange(-count, count+1)*dt
    values = np.r_[corr[size-count:], corr[:count+1]]
    k = int(np.argmax(values))
    far = np.abs(lag-lag[k]) > 1.0
    margin = float(values[k]-max(values[far])) if far.any() else 0.0
    return dict(ok=bool(k not in (0,len(values)-1) and values[k] > .25 and margin > .01),
                offset_s=float(lag[k]), peak=float(values[k]), margin=margin,
                boundary=bool(k in (0,len(values)-1)), reason='coarse_alignment_check')


def _diagnose_standard(ms, laser, settings, _edge=None):
    settings = dict(settings)
    original = dict(settings)
    allowed = {'tolerance_s', 'n_blocks', 'max_shift', 'peak_min', 'margin_min',
               'reject_tol', 'coarse_window_s', 'alignment_offset_s',
               'max_step_relative_deviation', 'min_coverage', 'closure_tolerance_ppm',
               'closure_iterations', 'scan_speed_um_s', 'line_change_limit_s',
               'exclude_record_boundary_blocks', 'skip_record_edge_lines'}
    if set(settings)-allowed:
        raise InputError('Unknown analysis setting(s): '+', '.join(sorted(set(settings)-allowed)))
    # The only non-numeric setting; it is validated here and restored below.
    exclude_boundary = settings.pop('exclude_record_boundary_blocks', False)
    if not isinstance(exclude_boundary, bool):
        raise InputError('exclude_record_boundary_blocks must be true or false')
    skip_edge = settings.pop('skip_record_edge_lines', True)
    if not isinstance(skip_edge, bool):
        raise InputError('skip_record_edge_lines must be true or false')
    defaults = dict(n_blocks=10, max_shift=25.0, peak_min=.25, margin_min=.01,
                    reject_tol=3.0, coarse_window_s=300., max_step_relative_deviation=.1,
                    min_coverage=.8, closure_tolerance_ppm=.2, closure_iterations=6,
                    line_change_limit_s=LINE_CHANGE_DEFAULT_LIMIT_S)
    settings = dict(defaults, **settings)
    if 'tolerance_s' not in settings or not math.isfinite(float(settings['tolerance_s'])) or float(settings['tolerance_s']) <= 0:
        raise InputError('An explicit positive tolerance_s is required for the intended use')
    for key, value in settings.items():
        if not isinstance(value, (int,float)) or isinstance(value,bool) or not math.isfinite(value):
            raise InputError('Analysis setting must be a finite number: '+key)
    for key in ('n_blocks', 'closure_iterations'):
        if int(settings[key]) != settings[key] or settings[key] < 1:
            raise InputError(key+' must be a positive integer')
        settings[key] = int(settings[key])
    if settings['n_blocks'] < 4 or settings['n_blocks'] > 100:
        raise InputError('n_blocks must be between 4 and 100')
    if settings['closure_iterations'] > 30:
        raise InputError('closure_iterations must not exceed 30')
    for key in ('max_shift','reject_tol','coarse_window_s','closure_tolerance_ppm','line_change_limit_s'):
        if settings[key] <= 0:
            raise InputError(key+' must be positive')
    settings['exclude_record_boundary_blocks'] = exclude_boundary
    settings['skip_record_edge_lines'] = skip_edge
    if not 0 < settings['min_coverage'] <= 1 or not 0 < settings['peak_min'] < 1 or not 0 < settings['margin_min'] < 1 or not 0 < settings['max_step_relative_deviation'] < 1:
        raise InputError('Invalid coverage, correlation gate or time-grid setting')
    if 'scan_speed_um_s' in settings and settings['scan_speed_um_s'] <= 0:
        raise InputError('scan_speed_um_s must be positive when supplied')
    t, y = np.asarray(ms.time, float), np.asarray(ms.signal,float)
    if t.ndim != 1 or y.shape != t.shape or len(t)<3 or not np.isfinite(t).all() or not np.isfinite(y).all() or not (np.diff(t)>0).all():
        raise InputError('Signal time and values must be finite 1D arrays with strictly increasing time')
    lines = laser.lines
    if not lines:
        raise InputError('No laser firing intervals')
    prev = -np.inf
    for row in lines:
        start, end = float(row['Start Time']), float(row['End Time'])
        if not (math.isfinite(start) and math.isfinite(end) and end > start and start >= prev):
            raise InputError('Laser firing intervals must be ordered, finite, non-overlapping and positive')
        prev = end
    dt = float(np.mean(np.diff(t)))
    grid_dev = float(np.max(abs(np.diff(t)/dt-1)))
    grid_error = float(np.max(abs(t-(t[0]+np.arange(len(t))*dt))))
    report = dict(status='insufficient_evidence', rate_ppm=None, settings=settings,
                  n_intervals=len(lines), n_signal_points=len(t), mean_spacing_s=dt,
                  max_step_relative_deviation=grid_dev, checks={},
                  index_grid_max_deviation_s=grid_error,
                  index_lag_error_bound_s=2*grid_error,
                  interpretation='Internal timing checks; not independent geometric accuracy or a clock attribution')
    # The laser-log screen needs only the logged times, so it is reported
    # whichever prerequisite or check stops the diagnosis below.
    screen = screen_line_changes(lines, settings['line_change_limit_s'])
    report['laser_log_screen'] = screen
    if grid_dev > settings['max_step_relative_deviation']:
        report['reason']='time_grid_not_supported_by_index_shift_estimator'
        return _native(report), []
    if 2*grid_error > settings['tolerance_s']/4:
        report['reason']='index_grid_approximation_exceeds_timing_budget'
        return _native(report), []
    if len(lines) < settings['n_blocks']*12:
        report['reason']='insufficient_laser_intervals'
        return _native(report), []
    if 'alignment_offset_s' in settings:
        alignment = dict(ok=True, offset_s=settings['alignment_offset_s'], source='explicit_user_offset')
    else:
        alignment = _coarse(t,y,lines,settings['coarse_window_s'])
        alignment['source']='automatic_candidate'
    report['alignment'] = alignment
    if not alignment['ok']:
        report['reason']='initial_alignment_not_identifiable'
        return _native(report), []
    aligned = _rescale(lines,0,alignment['offset_s'])
    # First pass: every line. A second pass (below) leaves out the lines at an end of the record
    # only when a block peak was found at a search boundary imposed by that end.
    edge = _edge or dict(start=0, end=0, applied=False, guard_s=None, sides=[])
    report['record_edge_lines'] = dict(enabled=skip_edge, applied=bool(edge['applied']), sides=list(edge.get('sides', [])),
                                       skipped_at_start=edge['start'], skipped_at_end=edge['end'], guard_s=edge['guard_s'],
                                       rule='When a block peak lies at a search boundary imposed by the start or end of the signal record, '
                                            'the diagnosis is repeated without the lines within the guard of that end of the record. '
                                            'The same lines are left out in every fit. Coverage is still measured against all lines.')
    whole = aligned
    aligned = whole[edge['start']:len(whole)-edge['end']]
    fit, df = _measure(t,y,aligned,settings)
    blocks = [] if df is None else df.to_dict(orient='records')
    retained = [] if fit is None else fit['retained']
    for b in blocks:
        reasons=[]
        if not b['peak'] > settings['peak_min']: reasons.append('low_peak_correlation')
        if not b['margin'] > settings['margin_min']: reasons.append('competing_peak')
        b['gate_pass']=not reasons
        b['retained']=b['block'] in retained
        if b.get('record_boundary_excluded'): reasons.append('record_boundary_excluded_by_option')
        elif not b['fit_eligible']: reasons.append('outside_signal_support')
        if b['gate_pass'] and b['fit_eligible'] and not b['retained'] and fit is not None: reasons.append('fit_outlier')
        b['reasons']=reasons
    report['computed_blocks']=len(blocks)
    report['gate_passing_blocks']=sum(b['gate_pass'] for b in blocks)
    if fit is None:
        report['reason']='fewer_than_four_usable_blocks'
        report['record_boundary_blocks']=_boundary_report(settings,df,None,None,None,[])
        return _native(report), _native(blocks)
    baseline_evidence=_evidence(fit,df,whole,settings)
    contradictory=baseline_evidence['contradictory_high_quality_blocks']
    for b in blocks:
        predicted=fit['intercept_s']+fit['rate_ppm']*1e-6*(b['t_mid']-fit['intercept_time_s'])
        b['model_residual_s']=b['lag_s']-predicted if b['fit_eligible'] else None
    report['contradictory_high_quality_blocks']=contradictory
    duration = float(lines[-1]['End Time']-lines[0]['Start Time'])
    coverage=baseline_evidence['coverage_fraction']
    boundary=baseline_evidence['peak_at_search_boundary']
    report['baseline']=fit
    report['duration_s']=duration
    report['coverage_fraction']=coverage
    trial=fit['rate_ppm']
    history=[]
    closed=None
    for iteration in range(settings['closure_iterations']):
        closed_whole = _rescale(whole,trial)
        closed_lines = closed_whole[edge['start']:len(closed_whole)-edge['end']]
        closed, closed_df = _measure(t,y,closed_lines,settings)
        residual=None if closed is None else closed['rate_ppm']
        history.append(dict(trial_ppm=trial,residual_ppm=residual))
        if residual is None or abs(residual)<settings['closure_tolerance_ppm']:
            break
        trial += residual
    converged=closed is not None and abs(closed['rate_ppm'])<settings['closure_tolerance_ppm']
    closure_evidence=_evidence(closed,closed_df,closed_whole,settings)
    closure_coverage=closure_evidence['coverage_fraction']
    closure_boundary=closure_evidence['peak_at_search_boundary']
    closure_contradictory=closure_evidence['contradictory_high_quality_blocks']
    report['closure']=dict(converged=converged, history=history,
                           residual_fit=closed,
                           residual_ppm=None if closed is None else closed['rate_ppm'],
                           coverage_fraction=closure_coverage,
                           peak_at_search_boundary=closure_boundary,
                           retained_blocks=[] if closed is None else closed['retained'],
                           contradictory_high_quality_blocks=closure_contradictory,
                           blocks=[] if closed_df is None else closed_df.to_dict(orient='records'))
    variants=[]
    for key, values in [('peak_min',(.20,.22,.28,.30)),('margin_min',(.015,.020)),
                        ('reject_tol',(2.,4.,5.)),('n_blocks',(8,12,15)),('max_shift',(20.,30.))]:
        for value in values:
            variant_settings=dict(settings,**{key:value})
            result,variant_df = _measure(t,y,aligned,settings,{key:value})
            evidence=_evidence(result,variant_df,whole,variant_settings)
            variants.append(dict(config=f'{key}={value}',rate_ppm=None if result is None else result['rate_ppm'],
                                 **evidence,blocks=[] if variant_df is None else variant_df.to_dict(orient='records')))
    for mode in ('all','trim'):
        result = _fit(df,settings['peak_min'],settings['margin_min'],settings['reject_tol'],mode)
        # Trimming intentionally removes the first and last qualified blocks.
        # Baseline/closure must cover the full record; this perturbation tests
        # the rate on an explicitly reduced domain, still retaining contradictions.
        evidence=_evidence(result,df,whole,settings,require_coverage=mode!='trim')
        variants.append(dict(config=mode,rate_ppm=None if result is None else result['rate_ppm'],
                             **evidence,blocks=df.to_dict(orient='records')))
    variants.append(dict(config='closure',rate_ppm=trial if converged else None,
                         **dict(closure_evidence,evidence_pass=converged and closure_evidence['evidence_pass']),
                         blocks=[] if closed_df is None else closed_df.to_dict(orient='records')))
    rates=[fit['rate_ppm']]+[v['rate_ppm'] for v in variants if v['rate_ppm'] is not None]
    envelope=(max(rates)-min(rates))*1e-6*duration
    tol=settings['tolerance_s']
    checks = dict(coverage=coverage>=settings['min_coverage'],
                  no_contradictory_high_quality_blocks=not contradictory,
                  peak_not_at_search_boundary=not boundary,
                  baseline_heldout=fit['heldout_max_s']<=tol,
                  all_sensitivity_configs_estimable=all(v['rate_ppm'] is not None for v in variants),
                  all_sensitivity_evidence_supported=all(v['evidence_pass'] for v in variants),
                  sensitivity_envelope=envelope<=tol,
                  closure_converged=converged,
                  closure_coverage=closure_coverage>=settings['min_coverage'],
                  closure_peak_not_at_boundary=not closure_boundary,
                  closure_no_contradictory_high_quality_blocks=not closure_contradictory,
                  closure_heldout=closed is not None and closed['heldout_max_s']<=tol,
                  closure_residual=closed is not None and abs(closed['rate_ppm'])*1e-6*duration<=tol,
                  no_flagged_line_change_interval=not screen['flagged'])

    boundary_report=_boundary_report(settings,df,fit,closed_df,closed,variants)
    excluded_blocks=boundary_report['excluded_blocks']
    # The option can also act only in an alternative setting of the stability
    # check (for example a finer block division); that still qualifies the result.
    option_acted=bool(excluded_blocks or boundary_report['excluded_in_sensitivity'])
    passed=all(checks.values())
    status='qc_failed' if not passed else ('supported_excluding_record_boundary_blocks' if option_acted else 'supported')
    report.update(checks=checks,sensitivity=variants,sensitivity_envelope_s=envelope,
                  rate_ppm=trial if converged else fit['rate_ppm'],
                  status=status,record_boundary_blocks=boundary_report,
                  failed_checks=[k for k,v in checks.items() if not v])
    if status=='supported_excluding_record_boundary_blocks':
        acted=[]
        if excluded_blocks:
            acted.append('block(s) '+', '.join(str(b) for b in excluded_blocks)+' left out of the reported fit')
        if boundary_report['excluded_in_sensitivity']:
            acted.append('a block left out in the stability alternatives '+', '.join(sorted(boundary_report['excluded_in_sensitivity'])))
        report['reason']=('All checks passed with the declared option exclude_record_boundary_blocks: '+'; '.join(acted)+'. '
                          'The search-edge criterion was not applied to an excluded block, and a lag beyond the record boundary cannot be seen.')
    report['predicted_end_to_end_lag_change_s']=report['rate_ppm']*1e-6*duration
    if 'scan_speed_um_s' in settings:
        report['predicted_end_to_end_displacement_um']=report['predicted_end_to_end_lag_change_s']*settings['scan_speed_um_s']
    if skip_edge and _edge is None:
        sides = _record_boundary_sides(df, closed_df, variants)
        if sides:
            candidate = record_edge_lines(t, whole, max(settings['max_shift'], EDGE_GUARD_MIN_S), settings['n_blocks'], tuple(sides))
            if candidate['applied']:
                candidate['sides'] = sorted(sides)
                return _diagnose_standard(ms, laser, original, _edge=candidate)
            report['record_edge_lines'].update(sides=sorted(sides), available_at_start=candidate['available_start'],
                                               available_at_end=candidate['available_end'],
                                               not_applied_reason='No line lies within the guard of that end, or fewer than twelve lines per block would remain')
    return _native(report), _native(blocks)


def diagnose(ms, laser, settings):
    from .short_acquisition import diagnose_short
    return diagnose_short(ms, laser, settings, _diagnose_standard)
