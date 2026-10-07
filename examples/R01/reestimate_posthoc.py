"""Re-estimate R01 using the documented post hoc boundary rule.

This example does not convert the rule into a default LA-LagCheck setting.
The signal choice Fe56 is fixed because only the authorised C13 and Fe56 values
are distributed; the earlier selection among all channels is not repeated.
"""
from pathlib import Path
import argparse, contextlib, io, json, sys
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from lagcheck import estimator as est
from lagcheck.models import SignalData, LaserData
import reconstruction as r
DROPPED=[]
def _measure_edge_drop(t, y, lines, settings, overrides=None):
    kw = {k: settings[k] for k in ('n_blocks', 'max_shift', 'peak_min', 'margin_min', 'reject_tol')}
    kw.update(overrides or {})
    with contextlib.redirect_stdout(io.StringIO()):
        _, df = est.core.analyse_loaded(t, y, lines, 'acquisition', **kw)
    if df is not None:
        dt = float((t[-1] - t[0]) / (len(t) - 1))
        pad = round(kw['max_shift'] / dt)
        edges = np.linspace(0, len(lines), kw['n_blocks'] + 1).astype(int)
        lower = []; upper = []; observed = []; clipped = []
        for row in df.to_dict(orient='records'):
            bi = int(row['block'])
            segment = lines[edges[bi]:edges[bi + 1]]
            i0 = int(np.searchsorted(t, segment[0]['Start Time']))
            i1 = int(np.searchsorted(t, segment[-1]['End Time']))
            j0 = max(0, i0 - pad); j1 = min(len(t), i1 + pad)
            lo = (j0 - i0) * dt; hi = (j1 - i1) * dt
            lower.append(lo); upper.append(hi)
            lag = row['lag_s']
            start = max(segment[0]['Start Time'], float(t[0]), float(t[0]) - lag)
            end = min(segment[-1]['End Time'], float(t[-1]), float(t[-1]) - lag)
            eligible = bool(np.isfinite(lag) and start <= segment[0]['Start Time'] + 1e-9
                            and end >= segment[-1]['End Time'] - 1e-9)
            cut_low = (i0 - pad) < 0 and np.isfinite(lag) and lag - lo <= 2 * dt
            cut_high = (i1 + pad) > len(t) and np.isfinite(lag) and hi - lag <= 2 * dt
            rec = bool(cut_low or cut_high)
            if rec:
                eligible = False; end = start          # Exclude this clipped block from fitting, coverage and agreement.
            clipped.append(dict(block=bi, side='start' if cut_low else ('end' if cut_high else ''), bound_s=lo if cut_low else hi) if rec else None)
            observed.append(dict(observed_start_s=start, observed_end_s=max(start, end),
                                 observed_n_lines=sum(l['Start Time'] >= start and l['End Time'] <= end for l in segment),
                                 fit_eligible=eligible))
        df['search_min_s'] = lower; df['search_max_s'] = upper
        df['near_search_boundary'] = (df.lag_s - df.search_min_s <= 2 * dt) | (df.search_max_s - df.lag_s <= 2 * dt)
        for key in ('observed_start_s', 'observed_end_s', 'observed_n_lines', 'fit_eligible'):
            df[key] = [row[key] for row in observed]
        df['record_clipped_search'] = [c is not None for c in clipped]
        fit = est._fit(df, kw['peak_min'], kw['margin_min'], kw['reject_tol'])
        if not overrides:
            for c, row in zip(clipped, df.to_dict(orient='records')):
                if c is None:
                    continue
                pred = None if fit is None else fit['intercept_s'] + fit['rate_ppm'] * 1e-6 * (row['t_mid'] - fit['intercept_time_s'])
                DROPPED.append(dict(block=c['block'], side=c['side'], bound_s=c['bound_s'], lag_s=row['lag_s'], predicted_s=pred,
                                    prediction_beyond_bound=None if pred is None else bool(pred <= c['bound_s'] if c['side'] == 'start' else pred >= c['bound_s'])))
        return fit, df
    return est._fit(df, kw['peak_min'], kw['margin_min'], kw['reject_tol']), df

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
out=args.output.resolve()
out.mkdir(parents=True,exist_ok=False)
finite=np.isfinite(r.y)&np.isfinite(r.t)
cuts=np.flatnonzero(np.diff(np.r_[False,finite,False].astype(int)))
segments=[(int(a),int(b)) for a,b in zip(cuts[::2],cuts[1::2]) if r.t[a]<r.L[-1]['End Time'] and r.t[b-1]>r.L[0]['Start Time']]
assert len(segments)==1
a,b=segments[0]
ms=SignalData(r.t[a:b],r.y[a:b],{'sample':'R01','finite_segment':[a,b],'excluded_signal_rows':len(r.t)-(b-a)})
laser=LaserData(r.L,{'sample':'R01'})
settings={'tolerance_s':.2}
standard,standard_blocks=est.diagnose(ms,laser,settings)
original_measure=est._measure
try:
    est._measure=_measure_edge_drop
    posthoc,posthoc_blocks=est.diagnose(ms,laser,settings)
finally:
    est._measure=original_measure
for name,report in [('preset_fe56',standard),('posthoc_fe56',posthoc)]:
    # Omit the synthetic arithmetic origin from all published result timestamps.
    def relative(x):
        if isinstance(x,dict):
            return {k:(v-r.BASE if k in {'t_mid','intercept_time_s','first_s','last_s','observed_start_s','observed_end_s'} and isinstance(v,(int,float)) else relative(v)) for k,v in x.items()}
        if isinstance(x,list): return [relative(v) for v in x]
        return x
    (out/(name+'.json')).write_text(json.dumps(relative(report),indent=2)+'\n',encoding='utf-8')
summary={'preset_fe56_status':standard['status'], 'posthoc_status':posthoc['status'],
         'posthoc_rate_ppm':posthoc.get('rate_ppm'), 'reference_rate_ppm':r.R_IT*1e6,
         'absolute_rate_difference_ppm':abs(posthoc['rate_ppm']-r.R_IT*1e6),
         'finite_segment':[a,b], 'excluded_signal_rows':len(r.t)-(b-a),
         'scope':'Post hoc Fe56 branch only; historical all-channel selection is not repeated. The final missing signal value is excluded from estimation but retained in the public input.'}
(out/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
assert summary['posthoc_status']=='supported'
assert summary['absolute_rate_difference_ppm']<1e-6
