"""Auditable visual inspection, independent of diagnostic acceptance."""
from pathlib import Path
from decimal import Decimal
import csv
import hashlib
import io
import json
import math
import textwrap
import numpy as np
from . import __version__
from .models import InputError, IntegrityError, SignalData
from .laser import load_laser
from .preview_geometry import _finite, _sha256, _read_json, _load_channel, _geometry, rasterize_true_y

MAX_CELLS = 1_000_000
MAX_AXIS = 4096


def _pyplot():
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        return plt
    except ImportError as exc:
        raise InputError('Default previews require matplotlib>=3.8; install the package requirements') from exc


def _map_lines(config, report, laser, ms):
    starts = np.array([v['Start Time'] for v in laser.lines])
    ends = np.array([v['End Time'] for v in laser.lines])
    fallback = report.get('alignment_fallback') or {}
    fit_report = fallback['diagnostic'] if fallback.get('adopted') else report
    alignment = fit_report.get('alignment') or {}
    offset, alignment_source = _finite(alignment.get('offset_s')), alignment.get('source')
    if fallback.get('adopted'):
        alignment_source = 'input_position_fallback; assumed input placement, NOT an identified clock origin'
    if offset is None:
        offset = _finite(config['analysis'].get('alignment_offset_s'))
        alignment_source = 'explicit_user_offset' if offset is not None else 'unregistered; zero display shift is not a known clock origin'
    display_offset = offset if offset is not None else 0.
    file = report.get('usable_file') or report.get('provisional_file')
    if file:
        provisional = not bool(report.get('usable_file'))
        receipt = report['provisional_export_verification' if provisional else 'export_verification']
        plan = report['provisional_export_plan' if provisional else 'export_plan']
        if _sha256(file) != receipt['output_sha256']:
            raise IntegrityError('Published log bytes changed before preview')
        written = load_laser(file, receipt['output_config'])
        if written.metadata['input_sha256'] != receipt['output_sha256']:
            raise IntegrityError('Published log changed while loading preview bytes')
        if len(written.lines) != len(laser.lines):
            raise InputError('Published log has a different firing-interval count')
        if plan['mode'] == 'affine':
            shift = float(written.table['origin']-Decimal(plan['coordinate_reference']['target_origin_seconds'])) + float(ms.time[0])
            alignment_source = 'published affine log in its explicitly declared MS coordinate'
        else:
            shift = display_offset + float(written.table['origin']-laser.table['origin'])
            alignment_source = 'published rate-only log with ' + str(alignment_source)
        mapped_starts = np.array([v['Start Time']+shift for v in written.lines])
        mapped_ends = np.array([v['End Time']+shift for v in written.lines])
        rate = float(receipt['rate_ppm'])
        source = 'published_provisional_log' if provisional else 'published_verified_log'
        provenance = dict(path=str(file), sha256=receipt['output_sha256'], mode=plan['mode'],
                          estimate_source='input_position_fallback' if fallback.get('adopted') else report.get('provisional_rate_source','acquisition_estimate'))
    else:
        rate = _finite(fit_report.get('rate_ppm'))
        source = ('input_position_fallback' if fallback.get('adopted') else 'acquisition_estimate') if rate is not None else None
        if rate is None:
            rate = _finite((config.get('export') or {}).get('reference_rate_ppm'))
            if rate is not None:
                source = 'user_reference_rate'
        provenance = None
        if rate is None:
            mapped_starts, mapped_ends = starts+display_offset, ends+display_offset
        else:
            scale = 1.+rate*1e-6
            if not math.isfinite(scale) or scale <= 0:
                raise InputError('Preview rate must have a finite positive time scale')
            constant = display_offset
            closure = fit_report.get('closure') or {}
            fit = closure.get('residual_fit')
            if closure.get('converged') and fit:
                mids = [b['t_mid'] for b in closure['blocks'] if b['block'] in fit['retained']]
                constant += fit['intercept_s'] + fit['rate_ppm']*1e-6*(float(np.median(mids))-fit['intercept_time_s'])
            elif source in ('acquisition_estimate','input_position_fallback') and fit_report.get('baseline'):
                fit = fit_report['baseline']
                constant += fit['intercept_s'] + rate*1e-6*(starts[0]+display_offset-fit['intercept_time_s'])
            mapped_starts = starts[0]+scale*(starts-starts[0])+constant
            mapped_ends = starts[0]+scale*(ends-starts[0])+constant
    lag = mapped_starts-starts
    return dict(rate=rate, source=source, provenance=provenance, alignment_offset_s=offset,
                alignment_source=alignment_source, starts=mapped_starts, ends=mapped_ends,
                lags=lag, constant=float(np.median(lag)) if rate is not None else display_offset,
                constant_source='median start-time shift of displayed mapping' if rate is not None else alignment_source)


def _grid(geometry):
    geometry = dict(geometry)
    x0 = float(min(geometry['along_start'].min(), geometry['along_end'].min()))
    x1 = float(max(geometry['along_start'].max(), geometry['along_end'].max()))
    y0, y1 = float(geometry['across'].min()), float(geometry['across'].max())
    pitch, row_pitch = geometry['pitch'], geometry['row_pitch']
    nx, ny = int(math.ceil((x1-x0)/pitch))+1, int(round((y1-y0)/row_pitch))+1
    factor = math.ceil(max(1., math.sqrt(nx*ny/MAX_CELLS), nx/MAX_AXIS, ny/MAX_AXIS))
    while True:
        px, py = pitch*factor, row_pitch*factor
        nx, ny = int(math.ceil((x1-x0)/px))+1, int(math.ceil((y1-y0)/py))+1
        if nx*ny <= MAX_CELLS and max(nx, ny) <= MAX_AXIS:
            break
        factor += 1
    geometry.update(pitch=px, row_pitch=py, downsample_factor=factor)
    return geometry, (x0, nx, y0, ny)


def _signal_view(png, ms, channel, headline, warnings):
    plt = _pyplot()
    t, y = np.asarray(ms.time), np.asarray(ms.signal)
    groups = np.array_split(np.arange(len(t)), min(2000, len(t)))
    selected=[]
    for group in groups:
        valid=group[np.isfinite(y[group])]
        if valid.size:
            selected.extend([valid[np.argmin(y[valid])],valid[np.argmax(y[valid])]])
    # Retain a NaN and its neighbouring points at every missing-data boundary,
    # so downsampling cannot connect across an omitted region.
    missing=~np.isfinite(y)
    edges=np.flatnonzero(np.diff(np.r_[False,missing,False].astype(int)))
    selected.extend(np.clip(np.r_[edges-1,edges],0,len(t)-1).tolist())
    indices=np.unique(selected)
    fig, ax = plt.subplots(figsize=(10, 4.8), layout='constrained')
    ax.plot(t[indices], y[indices], lw=.65)
    if missing.any():
        baseline=float(np.nanmin(y)) if np.isfinite(y).any() else 0.
        missing_times=t[np.flatnonzero(missing)[::max(1,int(missing.sum())//2000)]]
        ax.plot(missing_times,np.full(len(missing_times),baseline),'r|',markersize=9,label='missing/nonfinite signal (gap)')
        ax.legend(fontsize=8)
    ax.set(xlabel='MS time from recorded origin (s)', ylabel=channel,
           title='Raw signal inspection; no spatial reconstruction asserted')
    fig.suptitle(textwrap.fill(headline, 115), fontsize=10)
    fig.supxlabel('\n'.join(textwrap.fill(s, 125) for s in warnings), fontsize=8)
    fig.savefig(png, dpi=160, bbox_inches='tight', pad_inches=.15)
    plt.close(fig)
    return dict(view='signal_time', decimation='time-ordered per-bin minimum and maximum',
                samples_shown=len(indices), original_samples=len(t))


def render_partial_signal(run_dir, config, hashes, reason):
    """Inspect explicitly mapped numeric CSV values; never feed them to diagnosis."""
    cfg=config['ms']
    if cfg.get('format') != 'csv' or cfg.get('time_kind') != 'numeric':
        raise InputError('Partial-signal inspection requires an explicit numeric CSV mapping')
    units={'s':1.,'ms':.001,'us':1e-6,'min':60.}
    if cfg.get('time_unit') not in units:
        raise InputError('Partial-signal inspection needs a known explicit time unit')
    channels=cfg.get('signal_columns')
    if not isinstance(channels,list) or not channels or not cfg.get('time_column'):
        raise InputError('Partial-signal inspection needs explicit time and signal columns')
    path=Path(cfg['path']); raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != hashes['ms']:
        raise IntegrityError('MS bytes changed before partial-signal inspection')
    text=raw.decode(cfg.get('encoding','utf-8-sig'))
    start=cfg.get('header_row',0)
    if isinstance(start,bool) or not isinstance(start,int) or start<0:
        raise InputError('Partial-signal inspection requires a valid header_row')
    reader=csv.reader(io.StringIO(''.join(text.splitlines(True)[start:])),delimiter=cfg.get('delimiter',','),strict=True)
    header=next(reader)
    def column(name):
        found=[i for i,v in enumerate(header) if v.strip()==name.strip()]
        if len(found)!=1:
            raise InputError('Partial-signal inspection column must occur exactly once: '+str(name))
        return found[0]
    tc=column(cfg['time_column']); cols=[column(v) for v in channels]
    times=[]; values=[]; invalid=[]
    for number,row in enumerate(reader,1):
        try:
            stamp=float(row[tc])*units[cfg['time_unit']]
        except (IndexError,ValueError) as exc:
            raise InputError('Partial inspection cannot place a record with an invalid time coordinate') from exc
        if not math.isfinite(stamp):
            raise InputError('Partial inspection cannot place a nonfinite time coordinate')
        data=[]; bad=[]
        for name,col in zip(channels,cols):
            try: value=float(row[col])
            except (IndexError,ValueError): value=float('nan')
            if not math.isfinite(value):
                value=float('nan'); bad.append(name)
            data.append(value)
        if bad:
            invalid.append(dict(data_record=number,physical_line_end=start+reader.line_num,channels=bad,time_in_input_unit=row[tc]))
        times.append(stamp); values.append(data[0])
    time=np.asarray(times); signal=np.asarray(values)
    if len(time)<2 or not np.all(np.diff(time)>0) or not np.isfinite(signal).any():
        raise InputError('Partial-signal inspection needs ordered times and at least one finite signal value')
    ms=SignalData(time=time-time[0],signal=signal,metadata={})
    out=Path(run_dir); png=out/'preview_pair.png'; meta=out/'preview_pair.json'
    if png.exists() or meta.exists():
        raise InputError('Partial-signal preview already exists; refusing overwrite')
    warnings=['Strict input error: '+str(reason),
              f'{len(invalid)} source records contain missing or nonfinite requested signal values; no records were removed or repaired for diagnosis.',
              'Raw finite signal values only; missing regions remain gaps. No correction or spatial reconstruction.']
    for key,expected in hashes.items():
        if _sha256(config[key]['path']) != expected:
            raise IntegrityError('Input bytes changed before partial-signal inspection')
    view=_signal_view(png,ms,channels[0],'Strict input rejected; partial raw signal inspection',warnings)
    audit=dict(version=__version__,operation='inspection_only',analysis_cleaning_performed=False,
               input_sha256=hashes,source_path=str(path),invalid_signal_records=invalid,
               invalid_record_count=len(invalid),source_record_count=len(times),
               source_time_origin_s=float(time[0]),plot_time_basis='seconds from the first original MS record')
    with (out/'input_processing.json').open('x',encoding='utf-8') as f:
        json.dump(audit,f,ensure_ascii=False,indent=2,allow_nan=False)
    for key,expected in hashes.items():
        if _sha256(config[key]['path']) != expected:
            raise IntegrityError('Input bytes changed during partial-signal inspection')
    summary=dict(tool='la-lagcheck',version=__version__,png=str(png),json=str(meta),channel=channels[0],
                 rate_fitted=False,applied_rate_ppm=None,rate_source=None,report_status='input_or_tool_error',
                 inspection_only=True,input_processing=audit,warnings=warnings,**view)
    with meta.open('x',encoding='utf-8') as f:
        json.dump(summary,f,ensure_ascii=False,indent=2,allow_nan=False)
    return summary


def render_raw_diagnostic(run_dir, config, hashes, reason, name='preview_pair', details=None):
    out = Path(run_dir)
    png, meta = out/(name+'.png'), out/(name+'.json')
    if png.exists() or meta.exists():
        raise InputError('Preview already exists; refusing overwrite')
    channel, stored, note, ms = _load_channel(config, dict(input_sha256=hashes), None)
    for key, expected in hashes.items():
        if _sha256(config[key]['path']) != expected:
            raise IntegrityError('Input bytes changed; raw diagnostic preview refused')
    warnings = [str(reason), 'No correction or clock alignment is inferred from this raw signal view.']
    view = _signal_view(png, ms, channel, 'Input inspection', warnings)
    summary = dict(tool='la-lagcheck', version=__version__, png=str(png), json=str(meta),
                   channel=channel, rate_fitted=False, applied_rate_ppm=None, rate_source=None,
                   warnings=warnings, input_sha256=hashes, details=details, **view)
    for key, expected in hashes.items():
        if _sha256(config[key]['path']) != expected:
            raise IntegrityError('Input changed during raw diagnostic preview')
    with meta.open('x', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2, allow_nan=False)
    return summary


def render_preview_pair(run_dir, channel=None, name='preview_pair'):
    out = Path(run_dir).resolve()
    if not isinstance(name, str) or not name or Path(name).name != name or name.startswith('.'):
        raise InputError('name must be a plain file stem without directories')
    if (out/'INCOMPLETE.json').exists():
        raise InputError('This run is INCOMPLETE; only a separately marked raw diagnostic view is available')
    png, meta = out/(name+'.png'), out/(name+'.json')
    if png.exists() or meta.exists():
        raise InputError('Preview already exists; refusing overwrite')
    config, report = _read_json(out/'config.json'), _read_json(out/'report.json')
    try:
        if not all(_sha256(config[k]['path'])==v for k,v in report['input_sha256'].items()):
            raise IntegrityError('Input bytes changed before preview')
    except OSError as exc:
        raise IntegrityError('Input bytes are no longer available for preview') from exc
    channel, stored, channel_note, ms = _load_channel(config, report, channel)
    laser = load_laser(config['laser']['path'], {k:v for k,v in config['laser'].items() if k != 'path'})
    if laser.metadata['input_sha256'] != report['input_sha256']['laser']:
        raise IntegrityError('The laser log bytes changed before preview')
    mapping = _map_lines(config, report, laser, ms)
    geometry, speed, speed_source, uniform = _geometry(laser, ms, report.get('settings') or config['analysis'])
    if not geometry['source'].startswith('laser log'):
        durations = np.array([r['End Time']-r['Start Time'] for r in laser.lines])
        geometry.update(source='line order and time within firing intervals; no spatial coordinates',
                        unit='s', along_start=np.zeros(len(laser.lines)), along_end=durations,
                        across=np.arange(len(laser.lines), dtype=float), row_pitch=1.,
                        pitch=float(np.mean(np.diff(ms.time))), pitch_source='mean MS sample spacing')
    geometry, grid = _grid(geometry)
    x0, nx, y0, ny = grid
    lines = [dict(r, StartX=float(geometry['along_start'][i]), EndX=float(geometry['along_end'][i]),
                  StartY=float(geometry['across'][i])) for i, r in enumerate(laser.lines)]
    constant = mapping['constant']
    base = [dict(r, **{'Start Time': r['Start Time']+constant, 'End Time':r['End Time']+constant}) for r in lines]
    mapped = [dict(r, **{'Start Time':float(mapping['starts'][i]), 'End Time':float(mapping['ends'][i])}) for i,r in enumerate(lines)]
    unverified_peak = 'unverified_automatic_peak_candidate' in str(mapping['alignment_source'])
    base_title = ('unverified automatic alignment; constant only' if mapping['rate'] is not None else 'Unverified automatic alignment; no rate correction') if unverified_peak else ('one constant' if mapping['rate'] is not None else 'original / declared alignment')
    panel_lines = [(base_title, base)]
    if mapping['rate'] is not None:
        panel_lines.append(('correction: '+mapping['source'], mapped))
    images = [(title, rasterize_true_y(np.asarray(ms.time), np.asarray(ms.signal), rows, np.zeros(len(rows)),
                                      x0, nx, geometry['pitch'], y0, ny, geometry['row_pitch'])) for title,rows in panel_lines]
    warnings = []
    if report['status'] != 'supported':
        warnings.append('NOT VERIFIED: '+str(report.get('reason') or report['status']))
    failed = [k for k,v in report.get('checks', {}).items() if not v]
    failed += (report.get('post_export_verification') or {}).get('failed_checks', [])
    if failed:
        warnings.append('Failed checks: '+', '.join(dict.fromkeys(failed)))
    fallback = report.get('alignment_fallback') or {}
    if fallback.get('attempted'):
        warnings.append(fallback['note'])
        candidate_failed = [k for k,v in fallback['diagnostic'].get('checks',{}).items() if not v]
        if candidate_failed:
            warnings.append('Fallback candidate failed checks: '+', '.join(candidate_failed))
    warnings += laser.metadata.get('warnings', [])
    if unverified_peak:
        warnings.append('The automatic alignment candidate is not an identified clock origin.')
    if laser.metadata.get('format') == 'new_wave_laser_log':
        warnings.append('No verified firing speed: inherited positioning speeds are not used for spatial displacement estimates.')
    if geometry.get('note'):
        warnings.append(geometry['note'])
    if mapping['rate'] is None:
        warnings.append('No rate estimate or user reference rate: unverified automatic alignment only; no rate correction.' if unverified_peak else 'No rate estimate or user reference rate: only original/declared alignment shown.')
    if mapping['alignment_offset_s'] is None and (mapping['provenance'] or {}).get('mode') != 'affine':
        warnings.append('Clock origins are not aligned; zero is a display shift, not a verified origin.')
    if mapping['provenance'] is None and mapping['rate'] is not None:
        warnings.append('Estimated display mapping only; not an exported-log read-back.')
    if not report.get('correction_exported'):
        warnings.append('No verified corrected file is established by this view.')
    if geometry['downsample_factor'] > 1:
        warnings.append(f"Grid pitch enlarged {geometry['downsample_factor']} times; not full-resolution output.")
    post=report.get('post_export_verification')
    provisional_readback=report.get('provisional_readback')
    readback_status=('passed' if post['passed'] else 'failed') if post else ((provisional_readback or {}).get('status') or 'not_requested')
    export_status='verified' if report.get('correction_exported') else ('provisional_not_verified' if report.get('provisional_file') else 'not_published')
    stages=dict(input_diagnostic_status=report.get('input_diagnostic_status',report['status']),
                input_failed_checks=[k for k,v in report.get('checks',{}).items() if not v],
                readback_status=readback_status,readback_failed_checks=(post or provisional_readback or {}).get('failed_checks',[]),
                export_verification_status=export_status)
    headline = f"{laser.metadata.get('sample') or Path(config['laser']['path']).stem}; {channel}\nInput: {stages['input_diagnostic_status']}; read-back: {readback_status}; export: {export_status}"
    valid = [img[np.isfinite(img)] for _,img in images]
    if not any(v.size for v in valid):
        warnings.append('No MS samples fall within displayed firing intervals. Raw signal shown instead.')
        scale = _signal_view(png, ms, channel, headline, warnings)
        view = 'signal_time'
    else:
        view = 'stage_raster' if geometry['source'].startswith('laser log') else 'line_time'
        plt = _pyplot()
        values = np.log10(np.maximum(np.concatenate(valid), 1.))
        lo, hi = map(float, np.percentile(values, [2,98]))
        if hi-lo < .1:
            lo,hi = (hi+lo)/2-.05,(hi+lo)/2+.05
        footer='\n'.join(textwrap.fill(w,145) for w in warnings)
        header='\n'.join(textwrap.fill(line,115) for line in (headline+'\n'+str(mapping['alignment_source'])).splitlines())
        fig, axes = plt.subplots(1,len(images),figsize=(6*len(images),6+.16*(footer.count('\n')+header.count('\n'))),squeeze=False,layout='constrained')
        cmap=plt.get_cmap('viridis').copy(); cmap.set_bad('.85')
        from_log = view == 'stage_raster'
        extent=[x0,x0+nx*geometry['pitch'],y0+ny*geometry['row_pitch'],y0]
        for ax,(title,img) in zip(axes[0],images):
            with np.errstate(invalid='ignore'):
                data=np.log10(np.maximum(img,1.))
            artist=ax.imshow(data,cmap=cmap,vmin=lo,vmax=hi,extent=extent,aspect='equal' if from_log else 'auto',interpolation='nearest')
            ax.set(title=title,xlabel='recorded along-scan coordinate (um)' if from_log else 'time within firing interval (s)',
                   ylabel='recorded cross-scan coordinate (um)' if from_log else 'line order (not spatial distance)')
            if not np.isfinite(img).any():
                ax.text(.5,.5,'No samples in these intervals',transform=ax.transAxes,ha='center',bbox={'facecolor':'white'})
        fig.colorbar(artist,ax=axes[0].tolist(),shrink=.7,label='log10('+channel+'), floor 1')
        fig.suptitle(header,fontsize=10)
        fig.supxlabel(footer,fontsize=8)
        fig.savefig(png,dpi=160,bbox_inches='tight',pad_inches=.15); plt.close(fig)
        scale=dict(vmin=lo,vmax=hi,log10_percentiles=[2,98],extent=extent)
    new_wave_views = None
    if laser.metadata.get('new_wave'):
        from .new_wave_view import render_new_wave
        new_wave_views = render_new_wave(out, config, report, laser, ms, mapping, channel, _pyplot)
    if not all(_sha256(config[k]['path'])==v for k,v in report['input_sha256'].items()):
        raise IntegrityError('Input changed during preview; preview cannot be accepted')
    if mapping['provenance'] and _sha256(mapping['provenance']['path']) != mapping['provenance']['sha256']:
        raise IntegrityError('Published log changed during preview; preview cannot be accepted')
    starts=np.array([r['Start Time'] for r in laser.lines]); lag=mapping['lags']
    summary=dict(tool='la-lagcheck',version=__version__,png=str(png),json=str(meta),view=view,
                 report_status=report['status'],report_reason=report.get('reason'),warnings=warnings,
                 channel=channel,stored_channels=stored,channel_note=channel_note,
                 rate_fitted=_finite(report.get('rate_ppm')) is not None or bool(fallback.get('adopted')),applied_rate_ppm=mapping['rate'],
                 alignment_fallback=fallback,
                 rate_source=mapping['source'],mapping_provenance=mapping['provenance'],
                 baseline_rate_ppm=(report.get('baseline') or {}).get('rate_ppm'),report_rate_ppm=report.get('rate_ppm'),
                 alignment_offset_s=mapping['alignment_offset_s'],alignment_source=mapping['alignment_source'],
                 constant_lag_s=constant,constant_source=mapping['constant_source'],
                 first_line_lag_s=float(lag[0]) if mapping['rate'] is not None else None,
                 last_line_lag_s=float(lag[-1]) if mapping['rate'] is not None else None,
                 lag_span_s=float(lag[-1]-lag[0]) if mapping['rate'] is not None else None,
                 first_to_last_line_start_span_s=float(starts[-1]-starts[0]),
                 scan_speed_um_s=speed,scan_speed_source=speed_source,scan_speed_uniform_in_log=uniform,
                 geometry_source=geometry['source'],geometry_note=geometry['note'],axes_swapped=geometry['swapped'],
                 grid_shape=[ny,nx],pixel_pitch=geometry['pitch'],pixel_pitch_unit=geometry['unit'],row_pitch=geometry['row_pitch'],
                 grid_downsample_factor=geometry['downsample_factor'],colour_scale=scale,
                 panels=['raw signal timeline'] if view=='signal_time' else [t for t,_ in images],
                 mapping_displayed=view!='signal_time',
                 mapped_interval_starts_s=mapping['starts'].tolist(),mapped_interval_ends_s=mapping['ends'].tolist(),
                 input_processing=laser.metadata.get('input_processing'),input_sha256=report['input_sha256'],input_unchanged=True,
                 note='Visual inspection only; no independent spatial accuracy or receiving-importer validation.')
    if new_wave_views is not None:
        summary['new_wave_views'] = new_wave_views
    summary.update(stages)
    summary.update(original_input_check=report.get('original_input_check'),
                   processed_input_diagnostic_status=report.get('processed_input_diagnostic_status'))
    summary.update(first_line_displacement_um=float((lag[0]-constant)*speed) if speed is not None and mapping['rate'] is not None else None,
                   last_line_displacement_um=float((lag[-1]-constant)*speed) if speed is not None and mapping['rate'] is not None else None,
                   span_displacement_um=float((lag[-1]-lag[0])*speed) if speed is not None and mapping['rate'] is not None else None,
                   predicted_first_to_last_displacement_um=float(speed*mapping['rate']*1e-6*(starts[-1]-starts[0])) if speed is not None and mapping['rate'] is not None else None)
    with meta.open('x',encoding='utf-8') as f:
        json.dump(summary,f,ensure_ascii=False,indent=2,allow_nan=False)
    return summary
