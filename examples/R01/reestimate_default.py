"""Re-estimate R01 with the default LA-LagCheck settings.

The matching signal is chosen as the program does by default for CSV inputs:
the C13+Fe56 sum, C13 and Fe56 are compared by the peak of the automatic
starting-alignment search, and the candidate with the highest peak is used.
The diagnosis then runs with automatic starting alignment. The tolerance has
no default and is given here as 0.2 s. The rate, starting alignment,
iterated-fit intercept, retained blocks and block lags are compared with
reference_mapping.json and reference_blocks.csv, which define the mapping used
for Figure 3.

A second run declares the saved input placement (alignment_offset_s = 0) and
is reported for comparison. The earlier option exclude_record_boundary_blocks
still exists and is off by default; it is not needed for R01 because the
default rule skip_record_edge_lines leaves out the scan lines that start
within 30 s of the start of the signal record.
"""
from pathlib import Path
import argparse, json, sys
import numpy as np
import pandas as pd
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from lagcheck import estimator as est
from lagcheck.models import SignalData, LaserData
from lagcheck.signal_selection import select_signal

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
out = args.output.resolve()
out.mkdir(parents=True, exist_ok=False)

BASE = float(2**30)  # Synthetic arithmetic origin, not a real timestamp.
reference = json.loads((HERE/'reference_mapping.json').read_text(encoding='utf-8'))
reference_blocks = pd.read_csv(HERE/'reference_blocks.csv', float_precision='round_trip')
a = np.loadtxt(HERE/'R01_signals.csv', delimiter=',', skiprows=1)
t = a[:, 0] + BASE
columns = {'C13': a[:, 1], 'Fe56': a[:, 2]}
rows = pd.read_csv(HERE/'R01_workspace_intervals.csv', float_precision='round_trip').to_dict(orient='records')
lines = [dict(**{'Start Time': d['start_s'] + BASE, 'End Time': d['end_s'] + BASE},
              StartX=d['start_x_um'], EndX=d['end_x_um'], StartY=d['start_y_um'], EndY=d['end_y_um']) for d in rows]
finite = np.isfinite(t) & np.isfinite(columns['C13']) & np.isfinite(columns['Fe56'])
cuts = np.flatnonzero(np.diff(np.r_[False, finite, False].astype(int)))
segments = [(int(i), int(j)) for i, j in zip(cuts[::2], cuts[1::2])
            if t[i] < lines[-1]['End Time'] and t[j-1] > lines[0]['Start Time']]
assert len(segments) == 1
i0, i1 = segments[0]


def load(names):
    return SignalData(t[i0:i1], sum(columns[n][i0:i1] for n in names),
                      {'sample': 'R01', 'signal_columns': list(names), 'finite_segment': [i0, i1],
                       'excluded_signal_rows': len(t)-(i1-i0)})


TIME_KEYS = {'t_mid', 'intercept_time_s', 'first_s', 'last_s', 'observed_start_s', 'observed_end_s', 'next_start_time_s'}


def relative(x):
    # Omit the synthetic arithmetic origin from all published result timestamps.
    if isinstance(x, dict):
        return {k: (v-BASE if k in TIME_KEYS and isinstance(v, (int, float)) else relative(v)) for k, v in x.items()}
    if isinstance(x, list):
        return [relative(v) for v in x]
    return x


laser = LaserData(lines, {'sample': 'R01'})
results = {}
for name, extra in (('default_automatic', {}), ('default_declared_placement', {'alignment_offset_s': 0.0})):
    settings = dict({'tolerance_s': .2}, **extra)
    ms, selection = select_signal(load(['C13', 'Fe56']), laser, settings, load)
    report, blocks = est.diagnose(ms, laser, settings)
    report['signal_selection'] = selection
    (out/(name+'.json')).write_text(json.dumps(relative(report), indent=2)+'\n', encoding='utf-8')
    pd.DataFrame(relative(blocks)).to_csv(out/(name+'_blocks.csv'), index=False)
    results[name] = (report, blocks)

report, blocks = results['default_automatic']
declared = results['default_declared_placement'][0]
edge = report['record_edge_lines']
summary = {
    'status': report['status'], 'failed_checks': report['failed_checks'],
    'selected_signal': report['signal_selection']['selected_columns'],
    'candidate_alignment_peaks': {'+'.join(c['columns']): c['peak'] for c in report['signal_selection']['candidates']},
    'single_pass_rate_ppm': report['baseline']['rate_ppm'],
    'iterated_rate_ppm': report['rate_ppm'], 'reference_rate_ppm': reference['rate_ppm'],
    'absolute_rate_difference_ppm': abs(report['rate_ppm']-reference['rate_ppm']),
    'starting_alignment_s': report['alignment']['offset_s'],
    'iterated_fit_intercept_s': report['closure']['residual_fit']['intercept_s'],
    'scan_lines_left_out_at_start': edge['skipped_at_start'], 'scan_lines_left_out_at_end': edge['skipped_at_end'],
    'retained_blocks_after_iteration': report['closure']['retained_blocks'],
    'block_lags_s': [b['lag_s'] for b in blocks],
    'record_boundary_blocks': report['record_boundary_blocks']['baseline']+report['record_boundary_blocks']['closure'],
    'flagged_line_change_intervals': len(report['laser_log_screen']['flagged']),
    'declared_placement_status': declared['status'], 'declared_placement_failed_checks': declared['failed_checks'],
    'declared_placement_iterated_rate_ppm': declared['rate_ppm'],
    'finite_segment': [i0, i1], 'excluded_signal_rows': len(t)-(i1-i0),
    'scope': 'Default settings on the distributed C13 and Fe56 signals. The final missing signal value is excluded from '
             'estimation but retained in the public input. Internal timing checks; not independent spatial validation.'}
(out/'comparison.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
print(json.dumps(summary, indent=2))
assert summary['status'] == reference['status'] == 'supported' and not summary['failed_checks']
assert summary['selected_signal'] == [reference['channel']]
assert summary['absolute_rate_difference_ppm'] < 1e-6
assert abs(summary['starting_alignment_s']-reference['alignment']['offset_s']) < 1e-9
assert abs(summary['iterated_fit_intercept_s']-reference['closure']['residual_fit']['intercept_s']) < 1e-6
assert summary['scan_lines_left_out_at_start'] == reference['record_edge_lines']['skipped_at_start'] == 3
assert summary['scan_lines_left_out_at_end'] == 0
assert summary['retained_blocks_after_iteration'] == reference['closure']['retained_blocks'] == list(range(10))
assert np.allclose(summary['block_lags_s'], reference_blocks.lag_s, rtol=0, atol=1e-6)
assert np.allclose([b['t_mid']-BASE for b in blocks], reference_blocks.t_mid, rtol=0, atol=1e-5)
assert not summary['record_boundary_blocks']
assert summary['declared_placement_status'] == 'supported'
print('Default re-estimation matches reference_mapping.json and reference_blocks.csv.')
