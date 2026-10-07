"""Recorded-path/spot and temporal inspection for explicit New Wave inputs."""
import csv
import json
import math
import textwrap
from pathlib import Path
import numpy as np
from .models import InputError


def _write_csv(path, rows):
    if not rows:
        return
    with path.open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def render_new_wave(out, config, report, laser, ms, mapping, channel, pyplot):
    """Do not infer a tissue grid or connect different firing intervals."""
    out = Path(out)
    intervals = laser.metadata['new_wave']['intervals']
    t, signal = np.asarray(ms.time), np.asarray(ms.signal)
    constant = mapping['constant']
    unverified_peak = 'unverified_automatic_peak_candidate' in str(mapping['alignment_source'])
    first_title = ('Unverified automatic alignment; no rate correction' if mapping['rate'] is None else 'Unverified automatic alignment; constant only') if unverified_peak else 'original / constant placement'
    panels = [(first_title,
               np.array([v['Start Time'] + constant for v in laser.lines]),
               np.array([v['End Time'] + constant for v in laser.lines]))]
    if mapping['rate'] is not None:
        panels.append((('conditional mapping; origin unverified: ' if unverified_peak else 'conditional mapping: ') + str(mapping['source']), mapping['starts'], mapping['ends']))
    table, clouds, problems = [], [], []
    for panel, starts, ends in panels:
        cloud = []
        for item, start, end in zip(intervals, starts, ends):
            i0, i1 = np.searchsorted(t, [start, end])
            sample_times, values = t[i0:i1], signal[i0:i1]
            # New Wave repeated On rows may carry stale X/Y at the same time
            # as the first On. They establish state continuity, not a measured
            # intermediate path. Only firing endpoints define this display model.
            vertices = [item['event_records'][0], item['event_records'][-1]]
            try:
                vt = np.array([v['time_s'] for v in vertices], float)
                vx = np.array([float(v['fields']['X(um)']) for v in vertices])
                vy = np.array([float(v['fields']['Y(um)']) for v in vertices])
                if not np.isfinite(np.r_[vt, vx, vy]).all() or np.any(np.diff(vt) <= 0):
                    raise ValueError('nonfinite coordinates or non-increasing vertex times')
                stationary = np.ptp(vx) == 0 and np.ptp(vy) == 0
                mapped_vertex_times = start + (vt - item['start_s']) / (item['end_s'] - item['start_s']) * (end - start)
                xp = np.interp(sample_times, mapped_vertex_times, vx)
                yp = np.interp(sample_times, mapped_vertex_times, vy)
                if stationary:
                    cloud.append((item['sample'], [float(vx[0])], [float(vy[0])],
                                  [float(np.mean(values)) if len(values) else float('nan')], True))
                elif len(values):
                    cloud.append((item['sample'], xp, yp, values, False))
                geometry_error = None
            except (KeyError, ValueError) as exc:
                geometry_error = str(exc)
                stationary = None
                problems.append(dict(interval=item['interval'], panel=panel, reason=geometry_error))
                vx = vy = np.array([np.nan])
            table.append(dict(panel=panel, interval=item['interval'], sample=item['sample'], sequence_id=item['sequence_id'],
                source_start_record=item['start_record'], source_end_record=item['end_record'],
                mapped_start_s=float(start), mapped_end_s=float(end), sample_count=int(i1-i0),
                channel=channel, mean_signal=float(np.mean(values)) if len(values) else '',
                recorded_start_x_um=float(vx[0]) if np.isfinite(vx[0]) else '',
                recorded_start_y_um=float(vy[0]) if np.isfinite(vy[0]) else '',
                stationary=stationary, geometry_error=geometry_error or ''))
        clouds.append(cloud)
    all_stationary = all(c[-1] for cloud in clouds for c in cloud)
    labels = list(dict.fromkeys(i['sample'] for i in intervals))
    # Fixed spots are one position map. Moving acquisitions are separate panels
    # per exact sample label, preventing a distant standard from shrinking tissue.
    groups = ['all explicitly selected spots'] if all_stationary else labels
    plt = pyplot()
    fig, axes = plt.subplots(len(groups), len(panels), squeeze=False,
        figsize=(6*len(panels), 4.3*len(groups)+2), layout='constrained')
    valid = np.concatenate([np.asarray(v[3], float) for cloud in clouds for v in cloud]) if any(clouds) else np.array([])
    valid = valid[np.isfinite(valid)]
    lo, hi = np.percentile(np.log10(np.maximum(valid, 1)), [2, 98]) if valid.size else (0., 1.)
    if hi-lo < .1:
        lo, hi = lo-.05, hi+.05
    for row, group in enumerate(groups):
        for col, (panel, _, _) in enumerate(panels):
            ax = axes[row, col]
            artist = None
            for sample, x, y, values, stationary in clouds[col]:
                if not all_stationary and sample != group:
                    continue
                artist = ax.scatter(x, y, c=np.log10(np.maximum(values, 1)), vmin=lo, vmax=hi,
                    cmap='viridis', s=30 if stationary else 5, marker='o' if stationary else 's', linewidths=0)
            ax.set(title=group + '\n' + textwrap.fill(panel,55), xlabel='recorded X (um)', ylabel='recorded Y (um)')
            ax.set_aspect('equal', adjustable='datalim')
            if artist is not None:
                fig.colorbar(artist, ax=ax, shrink=.7, label='log10(' + channel + '), floor 1')
            else:
                ax.text(.5,.5,'No mapped finite observations',transform=ax.transAxes,ha='center')
    note = ('Spot colours are interval means; no spatial interpolation between spots.' if all_stationary else
            'Signal observations placed along the straight path between first-On and Off X/Y endpoints. No bridge between firing intervals or samples.')
    note += ' Intermediate On coordinates are state metadata, not verified trajectory vertices.'
    warning = 'NOT VERIFIED. ' + str(report.get('reason') or report['status'])
    alignment = 'Placement: ' + str(mapping['alignment_source'])
    fig.suptitle('New Wave recorded ' + ('spot positions' if all_stationary else 'scan paths') + '; ' + channel, fontsize=12)
    fig.supxlabel('\n'.join(textwrap.fill(v, 135) for v in [note, alignment, warning,
        'Spatial placement is a declared display model, not independent image accuracy or tested importer behaviour.']), fontsize=8)
    spatial = out/'new_wave_spatial.png'
    fig.savefig(spatial, dpi=160, bbox_inches='tight', pad_inches=.15)
    plt.close(fig)
    _write_csv(out/'new_wave_interval_signals.csv', table)

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), layout='constrained', sharex=True)
    # Min/max decimation retains short signal peaks in the whole-record display.
    indices = []
    for group in np.array_split(np.arange(len(t)), min(2500, len(t))):
        indices += [int(group[np.argmin(signal[group])]), int(group[np.argmax(signal[group])])]
    indices = np.unique(indices)
    axes[0].plot(t[indices], signal[indices], lw=.6)
    axes[0].set(ylabel=channel, title='Recorded signal and explicitly selected firing intervals')
    for n, (panel, starts, ends) in enumerate(panels):
        axes[1].broken_barh([(float(a), float(b-a)) for a,b in zip(starts,ends)], (n+.1,.65))
    axes[1].set(yticks=np.arange(len(panels))+.4, yticklabels=[p[0] for p in panels],
                xlabel='MS time from first recorded sample (s)')
    fig.supxlabel('\n'.join(textwrap.fill(v,145) for v in [alignment,warning]), fontsize=8)
    timeline=out/'new_wave_signal_time.png'
    fig.savefig(timeline,dpi=160,bbox_inches='tight',pad_inches=.15);plt.close(fig)

    blocks=json.loads((out/'blocks.json').read_text(encoding='utf-8'))
    fig,ax=plt.subplots(figsize=(10,5),layout='constrained')
    for kept,colour,label in [(True,'tab:blue','retained for fit'),(False,'tab:orange','not retained')]:
        rows=[b for b in blocks if b.get('retained')==kept and b.get('lag_s') is not None]
        ax.scatter([b['t_mid']/3600 for b in rows],[b['lag_s'] for b in rows],c=colour,label=label)
    fit=report.get('baseline')
    if fit:
        x=np.array([b['t_mid'] for b in blocks])
        ax.plot(x/3600,fit['intercept_s']+fit['rate_ppm']*1e-6*(x-fit['intercept_time_s']),lw=1)
    ax.set(xlabel='MS time (h)',ylabel='block lag (s)',title='Measured block lags; existing peak and margin gates retained')
    ax.legend(fontsize=8)
    fig.supxlabel(textwrap.fill(warning,120),fontsize=8)
    lagplot=out/'new_wave_block_lags.png'
    fig.savefig(lagplot,dpi=160,bbox_inches='tight',pad_inches=.15);plt.close(fig)
    summary=dict(adapter='new_wave_events_v1', interval_count=len(intervals), channel=channel,
        view='spot_position_map' if all_stationary else 'recorded_path_map',
        geometry_model=note, geometry_model_is_assumption=True, intermediate_on_positions_used=False,
        intermediate_on_explanation='Repeated On updates are preserved in report.json but can contain stale coordinates; no intermediate trajectory is asserted.',
        geometry_errors=problems, sample_groups=groups,
        source_sha256=report['input_sha256'], mapping_source=mapping['source'],
        mapping_provenance=mapping['provenance'], alignment_source=mapping['alignment_source'],
        input_diagnostic_status=report['status'], spatial_accuracy_verified=False,
        files=[str(p) for p in (spatial,timeline,lagplot,out/'new_wave_interval_signals.csv')],
        panels=[p[0] for p in panels], signal_values_modified=False)
    with (out/'new_wave_views.json').open('x',encoding='utf-8') as stream:
        json.dump(summary,stream,ensure_ascii=False,indent=2,allow_nan=False)
    return summary
