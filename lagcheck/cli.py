import argparse
import csv
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import sys
from . import __version__
from .models import InputError, IntegrityError
from .ms import load_signal
from .laser import load_laser, export_corrected, _time
from .estimator import diagnose
from .verification import export_plan, verify_written_log
from .fallback import try_alignment_fallback


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    with Path(path).open('x', encoding='utf-8') as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


def write_blocks(out, blocks, prefix='blocks'):
    write_json(out/(prefix+'.json'), blocks)
    if blocks:
        keys = sorted(set().union(*(r.keys() for r in blocks)))
        with (out/(prefix+'.csv')).open('x', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows({k: (json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v)
                         for k,v in row.items()} for row in blocks)


def _reject_known_corrected_input(path):
    input_hash = digest(path)
    for name in ('correction_manifest.json', 'provisional_manifest.json'):
        manifest = path.parent/name
        if not manifest.is_file():
            continue
        data = json.loads(manifest.read_text(encoding='utf-8'))
        if data.get('output_sha256') == input_hash:
            raise InputError('This file matches an earlier correction manifest (verified or provisional). Use the original laser log for export; '
                             'diagnosis without --export remains available. Do not apply an offset or rate twice.')


_EXPORT_KEYS = {'mode', 'ms_first_sample_at_laser_time', 'reference_note', 'provisional', 'reference_rate_ppm'}


def _validate_reference_rate(value):
    if value is None:
        return
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
            or 1 + value * 1e-6 <= 0):
        raise InputError('export.reference_rate_ppm must be a finite number in µs s⁻¹ giving a positive time-scale factor')


def _validate_export_options(options, laser):
    """Reject export configuration errors before any output is written (exit code 1)."""
    unknown = set(options) - _EXPORT_KEYS
    if unknown:
        raise InputError('Unknown export setting(s): ' + ', '.join(sorted(unknown)))
    mode = options.get('mode', 'rate_only')
    if mode not in ('rate_only', 'affine'):
        raise InputError('export.mode must be rate_only or affine')
    if not isinstance(options.get('provisional', True), bool):
        raise InputError('export.provisional must be true or false')
    _validate_reference_rate(options.get('reference_rate_ppm'))
    if mode == 'rate_only':
        if 'ms_first_sample_at_laser_time' in options:
            raise InputError('An explicit time reference requires export.mode=affine')
        return
    reference, note = options.get('ms_first_sample_at_laser_time'), options.get('reference_note')
    if not isinstance(reference, str) or not reference.strip() or not isinstance(note, str) or not note.strip():
        raise InputError('Affine export needs ms_first_sample_at_laser_time and a reference_note explaining its source; no clock origin is guessed')
    _, template = _time(reference, laser.metadata['config'])
    original_template = laser.table['templates'][0][0]
    if template is not None and ((template.tzinfo is None) != (original_template.tzinfo is None)):
        raise InputError('The time reference and laser timestamps must have consistent timezone awareness')


def _check_written_bytes(receipt, path):
    """The written file must still be the bytes the exporter checked."""
    if digest(path) != receipt['output_sha256']:
        raise InputError('Written file differs from the checked output bytes')


def _check_written_plan(receipt, rate, plan):
    """The written file must apply exactly the planned rate and constant."""
    if abs(float(receipt['rate_ppm']) - float(rate)) > 1e-9 * max(1.0, abs(float(rate))):
        raise InputError('Written rate differs from the export plan')
    if Decimal(str(receipt['offset_seconds'])) != Decimal(str(plan['offset_seconds'])):
        raise InputError('Written constant differs from the export plan')


def _check_written_mapping(receipt, rate, plan, path):
    """The written file must be the checked bytes and apply exactly the planned rate and constant."""
    _check_written_bytes(receipt, path)
    _check_written_plan(receipt, rate, plan)


def _markdown(report):
    lines = ['# Timing diagnostic report', '', f"Status: **{report['status']}**",
             f"Version: {__version__}", f"Rate estimate (µs s⁻¹): {report.get('rate_ppm')}",
             f"Verified corrected log available: {report['correction_exported']}",
             f"Usable file: {report.get('usable_file') or 'none'}",
             f"Provisional (not verified) log: {report.get('provisional_file') or 'none'}", '', report.get('reason', ''),
             '', '## Input diagnostic checks', '']
    lines.extend(f"- {k}: {'pass' if v else 'not met'}" for k,v in report.get('checks', {}).items())
    if report.get('short_acquisition'):
        short = report['short_acquisition']
        lines.extend(['', '## Conditional short-acquisition analysis', '', short['explanation'],
                      f"Blocks executed: {short['executed_blocks']}; usable fitted blocks: {short['usable_blocks']}",
                      'Unavailable sensitivity configurations: ' + ', '.join(short['unavailable_sensitivity_configs'])])
        lines.extend(f"- {k}: {v['state']}" + ('; ' + v['reason'] if v.get('reason') else '')
                     for k,v in short['check_execution'].items())
    fallback = report.get('alignment_fallback') or {}
    if fallback.get('attempted'):
        lines.extend(['', '## Separate input-placement fallback', '', fallback['note'],
                      'Candidate failed checks: ' + ', '.join(k for k,v in fallback['diagnostic'].get('checks', {}).items() if not v),
                      'Strict diagnostic fields above are unchanged.'])
    post = report.get('post_export_verification')
    if post:
        lines.extend(['', '## Actual exported CSV verification', '',
                      f"Residual estimate (µs s⁻¹): {post['residual_rate_ppm']}",
                      f"Single-pass residual (µs s⁻¹): {post['residual_baseline_rate_ppm']}",
                      f"Coverage fraction: {post['coverage_fraction']}",
                      f"Largest held-out residual (s): {post['heldout_max_s']}",
                      'The full diagnostic read the written CSV. No new automatic constant alignment was used.', ''])
        lines.extend(f"- {k}: {'pass' if v else 'not met'}" for k,v in post['checks'].items())
    if report.get('provisional_file'):
        readback = report.get('provisional_readback') or {}
        lines.extend(['', '## Provisional corrected log (not verified)', '',
                      f"File: {report['provisional_file']}",
                      f"Rate applied (µs s⁻¹): {report.get('provisional_rate_ppm')} ({report.get('provisional_rate_source')})",
                      f"Mapping: {report.get('provisional_mode')}",
                      f"Failed stage: {report.get('provisional_failure_stage')}",
                      'Failed checks: ' + (', '.join(report.get('provisional_failed_checks') or []) or str(report.get('reason') or 'see above')),
                      f"Fixed-mapping read-back: {readback.get('status')}; residual rate of the written file (µs s⁻¹): {readback.get('residual_rate_ppm')}",
                      (f"Read-back not run: {readback.get('reason')}" if readback.get('status') == 'not_run' else ''),
                      'The written timestamps were checked against the stated mapping; the timing evidence did not pass, '
                      'so this file is a provisional correction for the user to inspect, not a verified one.', ''])
    lines.extend(['', '## Reconstruction handoff', '', report.get('correction_note', 'No corrected file requested.'),
                  '', 'Original signal intensities and input files are unchanged.',
                  'Passing checks describes internal timing evidence, not independent spatial accuracy.',
                  'This command has not tested these bytes in iolite/HDIP. Use iolite-export for ESL formatting and timing checks, then verify the exact file in the receiving application.',
                  'Do not use candidate files from failed or INCOMPLETE runs.', ''])
    return '\n'.join(lines)


def _write_provisional(out, ms, laser, report, options, analysis, reference_rate_ppm=None):
    """Write a provisional (not verified) corrected log when the input checks did not all pass.

    The diagnosis is unchanged. The log applies the rate estimated from this
    acquisition, or a user-supplied reference rate when no rate could be
    estimated. The written timestamps are checked against the stated mapping,
    and the written file is re-diagnosed for information only. The returned
    manifest marks the file as provisional and lists the failed checks.
    """
    fallback = report.get('alignment_fallback') or {}
    fit_report = fallback['diagnostic'] if fallback.get('adopted') else report
    failed = [k for k, v in (fit_report.get('checks') or {}).items() if not v]
    rate, source = fit_report.get('rate_ppm'), 'input_position_fallback' if fallback.get('adopted') else 'acquisition_estimate'
    if rate is None and reference_rate_ppm is not None:
        rate, source = float(reference_rate_ppm), 'user_reference_rate'
    if rate is None:
        note = ('No provisional log: no rate could be estimated from this acquisition and no reference rate '
                'was given (--reference-rate-ppm). Resolve the failed condition, declare the starting alignment, '
                'or supply a reference rate measured on the same instrument pair.')
        report['provisional_note'] = note
        report['correction_note'] = 'No verified corrected log because the timing checks did not pass. ' + note
        return None
    requested = options.get('mode', 'rate_only')
    plan_note = ''
    closure = fit_report.get('closure') or {}
    affine_available = bool(closure.get('converged') and closure.get('residual_fit'))
    if source in ('acquisition_estimate', 'input_position_fallback') and (requested == 'rate_only' or affine_available):
        # Configuration errors were rejected earlier; any error here is not a reason to fall back.
        plan = export_plan(ms, laser, fit_report, options)
    else:
        if requested == 'affine':
            plan_note = (' An affine mapping needs the fitted constant of a converged iteration, which is not available here; '
                         'a rate-only provisional log was written instead.')
        plan = dict(mode='rate_only', offset_seconds='0', coordinate_reference=None,
                    initial_alignment_offset_s=(fit_report.get('alignment') or {}).get('offset_s'))
    candidate = out/'provisional_laser_log.csv'
    receipt = export_corrected(laser, candidate, rate, offset_seconds=plan['offset_seconds'])
    _check_written_mapping(receipt, rate, plan, candidate)
    readback = dict(status='not_run',
                    reason='No starting alignment is available for a fixed-mapping read-back; no new alignment is guessed.')
    if (fit_report.get('alignment') or {}).get('offset_s') is not None:
        # Byte or mapping errors raised here are execution errors, not informational read-back results.
        post, rerun, reblocks = verify_written_log(ms, laser, fit_report, candidate, receipt, plan, analysis)
        write_json(out/'provisional_readback_report.json', rerun)
        write_blocks(out, reblocks, 'provisional_readback_blocks')
        readback = dict(status='passed' if post['passed'] else 'timing_checks_failed',
                        residual_rate_ppm=post['residual_rate_ppm'],
                        residual_baseline_rate_ppm=post['residual_baseline_rate_ppm'],
                        checks=post['checks'], failed_checks=post['failed_checks'],
                        interpretation='Re-diagnosis of the written provisional file with the mapping held fixed')
    _check_written_mapping(receipt, rate, plan, candidate)
    what = ('the rate estimated at the original input placement' if fallback.get('adopted') else
            'the rate estimated from this acquisition' if source == 'acquisition_estimate' else 'the user-supplied reference rate')
    note = ('Provisional, NOT verified: the input timing checks did not all pass (' +
            (', '.join(failed) or str(report.get('reason') or 'see report')) + '). '
            f'This log applies {what} ({rate:.3f} µs s⁻¹) as a {plan["mode"]} mapping. '
            'Inspect the reconstruction before use and state the failed checks with any result.' + plan_note)
    if plan['mode'] == 'rate_only':
        note += ' Determine the constant offset in the reconstruction importer; do not apply the rate a second time.'
    if fallback.get('attempted'):
        note += ' ' + fallback['note']
    failure_stage = 'alignment_fallback' if fallback.get('adopted') else 'input_diagnostic'
    report.update(provisional_file=str(candidate), provisional_rate_ppm=rate, provisional_rate_source=source,
                  provisional_mode=plan['mode'], provisional_export_plan=plan,
                  provisional_export_verification=receipt, provisional_readback=readback,
                  provisional_failed_checks=failed, provisional_failure_stage=failure_stage,
                  provisional_note=note, correction_note=note)
    return dict(version=__version__, provisional=True, verified=False, mode=plan['mode'],
                rate_ppm=rate, rate_source=source, failed_checks=failed, failure_stage=failure_stage,
                alignment_fallback=fallback, input_diagnostic_status=report.get('input_diagnostic_status'),
                diagnostic_status=report['status'], reason=report.get('reason'), readback_status=readback['status'],
                output_file=str(candidate), output_sha256=receipt['output_sha256'], export_plan=plan,
                full_export_diagnostic_passed=False, importer_roundtrip_tested=False,
                instruction=note)


def _publish_provisional_from_candidate(out, report, candidate, receipt, plan, post):
    """After a failed read-back, publish the already checked candidate bytes as provisional.

    The correction is not applied again. A candidate whose written rate or
    constant differs from the plan is not published (the run keeps exit code
    2). Changed bytes are an execution error and are not downgraded: they
    raise, so the run ends with exit code 1 and INCOMPLETE.json.
    """
    _check_written_bytes(receipt, candidate)
    if digest(candidate) != post['diagnostic_input_sha256']:
        raise InputError('Candidate changed after the read-back diagnostic')
    try:
        _check_written_plan(receipt, report['rate_ppm'], plan)
    except InputError as exc:
        report['provisional_note'] = f'No provisional log: {exc}. The candidate is retained only for diagnosis.'
        return None
    payload = candidate.read_bytes()
    final = out/'provisional_laser_log.csv'
    with final.open('xb') as f:
        f.write(payload)
    if digest(final) != receipt['output_sha256']:
        raise InputError('Provisional file differs from the checked candidate bytes')
    failed = list(post['failed_checks'])
    note = ('Provisional, NOT verified: the input diagnostic passed, but the read-back of the written file did not ('
            + ', '.join(failed) + '). '
            f"This log keeps the checked candidate bytes ({float(report['rate_ppm']):.3f} µs s⁻¹, {plan['mode']} mapping); "
            'the correction is not applied again. Inspect the reconstruction before use and state the failed checks with any result.')
    if plan['mode'] == 'rate_only':
        note += ' Determine the constant offset in the reconstruction importer; do not apply the rate a second time.'
    report.update(provisional_file=str(final), provisional_rate_ppm=report['rate_ppm'],
                  provisional_rate_source='acquisition_estimate', provisional_mode=plan['mode'],
                  provisional_export_plan=plan, provisional_export_verification=dict(receipt, output_path=str(final)),
                  provisional_readback=dict(status='timing_checks_failed', residual_rate_ppm=post['residual_rate_ppm'],
                                            residual_baseline_rate_ppm=post['residual_baseline_rate_ppm'],
                                            checks=post['checks'], failed_checks=failed),
                  provisional_failed_checks=failed, provisional_failure_stage='readback',
                  provisional_note=note, correction_note=note)
    return dict(version=__version__, provisional=True, verified=False, mode=plan['mode'],
                rate_ppm=report['rate_ppm'], rate_source='acquisition_estimate', failed_checks=failed,
                failure_stage='readback', diagnostic_status=report['status'],
                input_diagnostic_status=report.get('input_diagnostic_status'), reason=report.get('reason'),
                readback_status='timing_checks_failed', output_file=str(final), output_sha256=receipt['output_sha256'],
                export_plan=plan, full_export_diagnostic_passed=False, importer_roundtrip_tested=False,
                instruction=note)


def _preview_after_run(out, channel):
    """Default picture drawn after report.json; a drawing failure is explicit."""
    try:
        from .preview import render_preview_pair
        summary = render_preview_pair(out, channel)
        return dict(status='written', png=summary['png'], json=summary['json'],
                    channel=summary['channel'], rate_fitted=summary['rate_fitted'])
    except Exception as exc:  # the diagnosis is complete; the preview must not undo it
        failure = dict(status='failed', error=f'{type(exc).__name__}: {exc}')
        if isinstance(exc,IntegrityError):
            write_json(out/'INCOMPLETE.json',dict(status='tool_or_input_error',error=str(exc),usable_file=None,
                       note='Byte integrity failed during preview; no output from this run is accepted'))
        else:
            try:
                from .preview import render_raw_diagnostic
                config=json.loads((out/'config.json').read_text(encoding='utf-8'))
                report=json.loads((out/'report.json').read_text(encoding='utf-8'))
                name='preview_raw' if (out/'preview_pair.png').exists() or (out/'preview_pair.json').exists() else 'preview_pair'
                raw=render_raw_diagnostic(out,config,report['input_sha256'],'Requested preview failed: '+str(exc),name=name)
                failure.update(status='fallback_written',png=raw['png'],json=raw['json'],channel=raw['channel'],rate_fitted=False)
            except Exception as fallback_exc:
                failure['fallback_error']=str(fallback_exc)
                if isinstance(fallback_exc,IntegrityError) and not (out/'INCOMPLETE.json').exists():
                    write_json(out/'INCOMPLETE.json',dict(status='tool_or_input_error',error=str(fallback_exc),usable_file=None))
        write_json(out/'preview_error.json', failure)
        print(json.dumps(dict(preview=failure), ensure_ascii=False), file=sys.stderr)
        return failure


def _raw_error_preview(out, config, before, exc, partial_signal=False):
    """Keep the execution error and inspect original bytes only, never candidate outputs."""
    try:
        from .preview import render_raw_diagnostic, render_partial_signal
        if any(digest(config[k]['path']) != v for k,v in before.items()):
            return
        out.mkdir(parents=True, exist_ok=True)
        if not (out/'config.json').exists():
            write_json(out/'config.json', config)
        if not (out/'INCOMPLETE.json').exists():
            write_json(out/'INCOMPLETE.json', dict(status='tool_or_input_error', error=str(exc), usable_file=None))
        render=render_partial_signal if partial_signal else render_raw_diagnostic
        render(out, config, before, 'Execution error: '+str(exc))
    except Exception as preview_exc:
        print(json.dumps(dict(preview=dict(status='failed', error=str(preview_exc))),ensure_ascii=False),file=sys.stderr)


def _run_segments(config, out, before, segmentation, export, preview, options, alignment_fallback):
    from .preview import render_raw_diagnostic
    out.mkdir(parents=True, exist_ok=False)
    write_json(out/'config.json',config)
    entries=[]
    for item in segmentation['segments']:
        dest=out/f"segment_{item['segment']:03d}"
        child=json.loads(json.dumps(config))
        child['laser']['record_range']=[item['start_record'], item['end_record']]
        if item['sample']:
            child['laser']['sample_selector']=item['sample']
        if 'alignment_offset_s' in child['analysis']:
            child['analysis']['alignment_offset_s'] += item['origin_shift_s']
        entry=dict(item,output=str(dest),declared_alignment_offset_s=child['analysis'].get('alignment_offset_s'))
        if item.get('on_record_count') == 0 and child['laser'].get('format','laser_state_csv') == 'laser_state_csv':
            dest.mkdir()
            note='No firing On records in this fragment; metadata retained, no tissue image or corrected log is invented.'
            fragment=dict(status='no_firing_fragment',reason=note,source=item)
            write_json(dest/'report.json',fragment)
            entry.update(exit_code=2,status='no_firing_fragment',reason=note,preview=None)
        else:
            try:
                code, report=run(child,dest,export,preview=preview,alignment_fallback=alignment_fallback,
                                 _fallback_input_offset_s=item['origin_shift_s'])
                entry.update(exit_code=code,status=report['status'],rate_ppm=report.get('rate_ppm'),
                             preview=report.get('preview'),provisional_file=report.get('provisional_file'),usable_file=report.get('usable_file'),
                             alignment_fallback=report.get('alignment_fallback'))
            except (InputError,ValueError,KeyError,OSError) as exc:
                entry.update(exit_code=1,status='input_or_tool_error',error=str(exc))
        entries.append(entry)
    if any(digest(config[k]['path']) != v for k,v in before.items()):
        write_json(out/'INCOMPLETE.json',dict(status='tool_or_input_error',error='Input changed during segmented run'))
        raise InputError('Input changed during segmented run')
    code=1 if any(v['exit_code']==1 for v in entries) else 2
    report=dict(tool='la-lagcheck',version=__version__,status='segmented_with_errors' if code==1 else 'segmented_not_jointly_verified',
                reason='Contiguous sample/sequence fragments were processed separately; no joint correction was made.',
                rate_ppm=None,correction_exported=False,usable_file=None,provisional_file=None,
                exported_log_diagnostic_rerun=False,input_sha256=before,segments=entries)
    original_check=dict(status='requires_processing',reasons=['sample_or_sequence_segmentation'],
                        assessment='Source-record inspection; fragments are not one joint acquisition')
    report['original_input_check']=original_check
    segmentation['original_input_check']=original_check
    write_json(out/'segments.json',dict(segmentation,segments=entries))
    write_json(out/'input_processing.json',dict(segmentation,segments=entries))
    write_json(out/'report.json',report)
    (out/'report.md').write_text('# Segmented input\n\n'+report['reason']+'\n\n'+
        '\n'.join(f"- Segment {v['segment']}, records {v['start_record']} to {v['end_record']}: {v['status']}" for v in entries),encoding='utf-8')
    if preview is not False:
        try:
            fragment_summary='; '.join(f"segment {v['segment']}: {v['status']}" for v in entries)
            summary=render_raw_diagnostic(out,config,before,report['reason']+' '+fragment_summary+'. See segments.json and per-fragment previews.',details=entries)
            report['preview']=dict(status='written',png=summary['png'],json=summary['json'],channel=summary['channel'],rate_fitted=False)
        except Exception as exc:
            report['preview']=dict(status='failed',error=str(exc)); code=1
            write_json(out/'preview_error.json',report['preview'])
    return code,report


def run(config_path, output, export=False, export_mode=None, ms_first_sample_time=None, reference_note=None, preview=None,
        provisional=None, reference_rate_ppm=None, alignment_fallback=True, _fallback_input_offset_s=0.):
    if isinstance(config_path, dict):
        config = json.loads(json.dumps(config_path))
        base = Path.cwd()
    else:
        config_path = Path(config_path).resolve()
        config = json.loads(config_path.read_text(encoding='utf-8-sig'))
        base = config_path.parent
    if set(config)-{'ms','laser','analysis','export'} or any(k not in config for k in ('ms','laser','analysis')):
        raise InputError('Configuration requires ms, laser and analysis, with optional export settings')
    options = dict(config.get('export', {}))
    for key, value in [('mode', export_mode), ('ms_first_sample_at_laser_time', ms_first_sample_time), ('reference_note', reference_note),
                       ('provisional', provisional), ('reference_rate_ppm', reference_rate_ppm)]:
        if value is not None:
            options[key] = value
    config['export'] = options
    config['laser'].setdefault('zero_duration_policy','exclude')
    # Provisional output (new in 0.2.0) is on by default with --export; it never
    # changes the diagnosis or the verified export below.
    provisional = options.get('provisional', True)
    reference_rate_ppm = options.get('reference_rate_ppm')
    if not isinstance(provisional, bool):
        raise InputError('export.provisional must be true or false')
    _validate_reference_rate(reference_rate_ppm)
    plan_options = {k: v for k, v in options.items() if k not in ('provisional', 'reference_rate_ppm')}
    paths = {}
    for kind in ('ms','laser'):
        if not isinstance(config[kind], dict) or not isinstance(config[kind].get('path'), str):
            raise InputError(kind+'.path must be a filename')
        paths[kind] = (base/config[kind]['path']).resolve()
        if not paths[kind].is_file():
            if kind=='laser' and preview is not False and paths.get('ms') and paths['ms'].is_file():
                raw_out=Path(output).resolve()
                if not raw_out.exists():
                    _raw_error_preview(raw_out,config,{'ms':digest(paths['ms'])},InputError('laser source does not exist'))
            raise InputError(kind+' source does not exist')
        config[kind]['path'] = str(paths[kind])
    out = Path(output).resolve()
    if out.exists():
        raise InputError('Output directory must not already exist; originals/results are never overwritten')
    if export:
        _reject_known_corrected_input(paths['laser'])
    before = {k: digest(v) for k,v in paths.items()}
    ms_config = {k:v for k,v in config['ms'].items() if k!='path'}
    try:
        if ms_config.get('format')=='netcdf':
            from .netcdf import load_netcdf_signal
            ms = load_netcdf_signal(paths['ms'], ms_config)
        else:
            ms = load_signal(paths['ms'], ms_config)
    except (InputError,ValueError,OSError) as exc:
        if preview is not False:
            _raw_error_preview(out,config,before,exc,partial_signal=True)
        raise
    try:
        from .processing import inspect_segments
        segmentation=inspect_segments(paths['laser'],{k:v for k,v in config['laser'].items() if k!='path'})
        if segmentation['source_sha256'] != before['laser']:
            raise IntegrityError('Laser bytes read for segmentation differ from the initial input snapshot')
        if len(segmentation['segments']) > 1:
            return _run_segments(config,out,before,segmentation,export,preview,options,alignment_fallback)
        laser = load_laser(paths['laser'], {k:v for k,v in config['laser'].items() if k!='path'})
    except (InputError,ValueError,KeyError,OSError) as exc:
        if preview is not False:
            _raw_error_preview(out,config,before,exc)
        raise
    config['laser'].update(laser.metadata['config'])
    if export:
        _validate_export_options(options, laser)
    try:
        report, blocks = diagnose(ms, laser, config['analysis'])
        fallback, fallback_blocks = try_alignment_fallback(ms, laser, report, config['analysis'], diagnose,
            enabled=alignment_fallback, provisional=provisional, input_offset_s=_fallback_input_offset_s)
    except (InputError,ValueError,KeyError) as exc:
        if preview is not False:
            _raw_error_preview(out,config,before,exc)
        raise
    processed_status=report['status']
    processing=laser.metadata['input_processing']
    required=[]
    if processing['excluded_count']:
        required.append('zero_duration_intervals')
    if config['laser'].get('record_range') is not None:
        required.append('sample_or_sequence_segmentation')
    original_check=dict(status='requires_processing' if required else 'accepted',reasons=required,
                        assessment='Source-record inspection; not a second run of the timing diagnostic')
    processing.update(original_input_check=original_check,processed_input_diagnostic_status=processed_status,
                      publication_policy='processed inputs are provisional only' if required else 'normal two-stage verification')
    if required and processed_status=='supported':
        report['status']='processed_input_diagnostic_supported'
        report['reason']='The processed intervals support the model; the original input required '+', '.join(required)
    report.update(original_input_check=original_check,processed_input_diagnostic_status=processed_status)
    report['alignment_fallback'] = fallback
    report.update(tool='la-lagcheck', version=__version__, input_sha256=before,
                  input_metadata={'ms':ms.metadata, 'laser':laser.metadata},
                  input_diagnostic_status=report['status'], correction_exported=False,
                  importer_roundtrip_tested=False, exported_log_diagnostic_rerun=False,
                  receiver_interval_trajectory_agreement_verified=False,
                  usable_file=None, candidate_file=None)
    if any(digest(paths[k])!=v for k,v in before.items()):
        raise InputError('Input changed during analysis; no correction will be exported')
    out.mkdir(parents=True, exist_ok=False)
    provisional_manifest = None
    try:
        write_json(out/'config.json', config)
        write_json(out/'input_processing.json', laser.metadata['input_processing'])
        write_blocks(out, blocks)
        if fallback['attempted']:
            write_json(out/'alignment_fallback_report.json', fallback['diagnostic'])
            write_blocks(out, fallback_blocks, 'alignment_fallback_blocks')
        if export and report['status']=='supported':
            plan = export_plan(ms, laser, report, plan_options)
            report['export_plan'] = plan
            candidate = out/'candidate_laser_log.csv'
            receipt = export_corrected(laser, candidate, report['rate_ppm'], offset_seconds=plan['offset_seconds'])
            report['candidate_file'] = str(candidate)
            report['export_verification'] = receipt
            post, rerun, reblocks = verify_written_log(ms, laser, report, candidate, receipt, plan, config['analysis'])
            report['exported_log_diagnostic_rerun'] = True
            report['post_export_verification'] = post
            write_json(out/'post_export_report.json', rerun)
            write_blocks(out, reblocks, 'post_export_blocks')
            if post['passed']:
                # Publish exactly the fully diagnosed on-disk bytes. Retain the
                # candidate for the audit; no file is moved or overwritten.
                payload = candidate.read_bytes()
                if hashlib.sha256(payload).hexdigest()!=post['diagnostic_input_sha256']:
                    raise InputError('Candidate changed before publication')
                final = out/'corrected_laser_log.csv'
                with final.open('xb') as f:
                    f.write(payload)
                final_parsed = load_laser(final, receipt['output_config'])
                if digest(final)!=post['diagnostic_input_sha256'] or final_parsed.metadata['input_sha256']!=post['diagnostic_input_sha256']:
                    raise InputError('Published CSV differs from the fully diagnosed bytes')
                receipt.update(output_path=str(final), output_sha256=digest(final))
                post['published_identical_bytes'] = True
                report['correction_exported'] = True
                report['usable_file'] = str(final)
                report['correction_note'] = (plan.get('import_instruction') if plan['mode']=='affine' else
                    'Rate-only corrected log: determine its constant offset in the reconstruction importer. '
                    'Before image use, verify that firing-interval endpoints and continuous-trajectory transition times '
                    'use the same corrected time coordinate after the offset. '
                    'Do not apply the original rate a second time. Importer compatibility remains untested.')
            else:
                report['status'] = 'export_qc_failed'
                report['reason'] = 'Actual exported CSV did not pass all residual timing checks'
                report['correction_note'] = 'No usable corrected file. The candidate CSV is retained only for diagnosis; do not import it as an accepted correction.'
                if provisional:
                    provisional_manifest = _publish_provisional_from_candidate(out, report, candidate, receipt, plan, post)
        elif export and provisional:
            provisional_manifest = _write_provisional(out, ms, laser, report, plan_options, config['analysis'], reference_rate_ppm)
        elif export:
            report['correction_note'] = 'No corrected log exported because timing checks did not pass'
        report['input_unchanged'] = all(digest(paths[k])==v for k,v in before.items())
        if not report['input_unchanged']:
            raise InputError('Input changed before completion')
        write_json(out/'report.json', report)
        (out/'report.md').write_text(_markdown(report), encoding='utf-8')
        if provisional_manifest is not None:
            # Written last, like the verified manifest; it marks the file as provisional.
            write_json(out/'provisional_manifest.json', dict(provisional_manifest, input_sha256=before,
                       acquisition_identity=laser.metadata['acquisition_identity']))
        if report['correction_exported']:
            # The manifest is the last publication marker, after input hashes
            # and both human/machine-readable final reports are complete.
            write_json(out/'correction_manifest.json', dict(version=__version__, mode=plan['mode'],
                       input_sha256=before, output_sha256=digest(final), output_file=str(final),
                       export_plan=plan, full_export_diagnostic_passed=True,
                       acquisition_identity=laser.metadata['acquisition_identity'],
                       receiver_interval_trajectory_agreement_verified=False,
                       importer_roundtrip_tested=False, instruction=report['correction_note']))
        code = 0 if report['status']=='supported' else 2
    except Exception as exc:
        write_json(out/'INCOMPLETE.json', dict(status='tool_or_input_error', error=str(exc),
                   usable_file=None, note='Do not use any outputs from this incomplete run'))
        if preview is not False:
            _raw_error_preview(out,config,before,exc)
        raise
    if preview is not False:
        # Drawn only after report.json and any manifest exist. It reads the run
        # directory back; the written report, blocks and logs are not touched.
        report['preview'] = _preview_after_run(out, None if preview in (True,None) else preview)
        if report['preview']['status'] != 'written':
            code=1
    return code, report


def main(argv=None):
    parser = argparse.ArgumentParser(description='LA-LagCheck: per-acquisition laser/MS lag rate, checks and verified corrected log')
    parser.add_argument('--version', action='version', version=__version__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('run', help='Diagnose and optionally correct one acquisition from a JSON configuration')
    p.add_argument('config'); p.add_argument('--output', required=True)
    p.add_argument('--export', action='store_true', help='Export and fully re-diagnose the written CSV; publish only after checks pass')
    p.add_argument('--export-mode', choices=('rate_only','affine'))
    p.add_argument('--ms-first-sample-time', help='Explicit time of first stored MS sample in the laser output timestamp coordinate (affine only)')
    p.add_argument('--reference-note', help='Where that first-sample time reference was verified')
    p.add_argument('--no-provisional', dest='provisional', action='store_false', default=None,
                   help='Disable alignment fallback and provisional (not verified) log export')
    p.add_argument('--no-alignment-fallback', dest='alignment_fallback', action='store_false', default=True,
                   help='Do not estimate a separate provisional rate at the input placement after eligible automatic-alignment failures')
    p.add_argument('--reference-rate-ppm', type=float, default=None,
                   help='Rate in µs s⁻¹ to use for a provisional log when no rate can be estimated from this acquisition '
                        '(for example the median rate measured on the same instrument pair)')
    p.add_argument('--preview', nargs='?', const=True, default=None, metavar='CHANNEL',
                   help='Select the channel for the default preview_pair.png. Actual published mappings take precedence; '
                        'otherwise show the current estimated or explicitly supplied rate, or raw alignment when no rate exists. '
                        'The first stored channel is the default. Failed checks remain visible; --no-preview disables drawing.')
    p.add_argument('--no-preview', dest='preview', action='store_const', const=False,
                   help='Explicitly disable the default visual preview')
    p = sub.add_parser('wizard', help='Guided MassHunter/laser-file selection without writing JSON')
    p.add_argument('--no-alignment-fallback',dest='alignment_fallback',action='store_false',default=True)
    p.add_argument('--no-provisional',dest='provisional',action='store_false',default=None)
    p.add_argument('--no-preview',dest='preview',action='store_false',default=True,
                   help='Disable the default visual preview')
    p = sub.add_parser('preview', help='Draw the preview pair for a completed run directory without changing its results')
    p.add_argument('run_dir', help='Directory holding config.json and report.json of a completed run')
    p.add_argument('--channel', help='Signal column or variable to draw; default: the first stored channel')
    p.add_argument('--name', default='preview_pair', help='File stem inside the run directory (default preview_pair); existing files are never overwritten')
    p = sub.add_parser('iolite-export', help='Convert a completed canonical run to ESL millisecond CSV, then recheck saved bytes under its fixed timing mapping')
    p.add_argument('canonical_run', help='Directory containing config.json, report.json, correction_manifest.json and corrected_laser_log.csv')
    p.add_argument('--output', required=True, help='New result directory; existing directories are rejected')
    p = sub.add_parser('inspect-netcdf', help='List classic NetCDF metadata without guessing a mapping')
    p.add_argument('path')
    args = parser.parse_args(argv)
    try:
        if args.command=='iolite-export':
            from .iolite import export_iolite
            code, report = export_iolite(args.canonical_run, args.output)
            print(json.dumps({k: report[k] for k in ('status', 'version', 'canonical_mode',
                             'correction_exported', 'usable_file', 'provisional', 'provisional_file',
                             'importer_roundtrip_tested', 'gui_workflow_tested')}, ensure_ascii=False))
            return code
        if args.command=='inspect-netcdf':
            from .netcdf import inspect_netcdf
            print(json.dumps(inspect_netcdf(args.path), ensure_ascii=False, indent=2, allow_nan=False))
            return 0
        if args.command=='preview':
            from .preview import render_preview_pair
            summary = render_preview_pair(args.run_dir, args.channel, args.name)
            print(json.dumps({k: summary[k] for k in ('channel', 'rate_fitted', 'constant_lag_s',
                             'baseline_rate_ppm', 'png', 'json')}, ensure_ascii=False))
            return 0
        if args.command=='wizard':
            from .wizard import collect
            config, output = collect()
            code, report = run(config, output, True, preview=args.preview,provisional=args.provisional,
                               alignment_fallback=args.alignment_fallback)
        else:
            output = args.output
            code, report = run(args.config, output, args.export, args.export_mode,
                               args.ms_first_sample_time, args.reference_note, args.preview,
                               args.provisional, args.reference_rate_ppm, args.alignment_fallback)
        print(json.dumps({'status':report['status'], 'rate_ppm':report.get('rate_ppm'),
                          'rate_unit':'µs s⁻¹', 'correction_exported':report['correction_exported'],
                          'exported_log_diagnostic_rerun':report['exported_log_diagnostic_rerun'],
                          'usable_file':report['usable_file'], 'provisional_file':report.get('provisional_file'),
                          'output':str(Path(output).resolve()),
                          'preview':report.get('preview')}, ensure_ascii=False))
        if args.command=='wizard':
            if report['correction_exported']:
                print('Verified log: '+str(report['usable_file']))
            elif report.get('provisional_file'):
                print('Provisional (not verified) log: '+str(report['provisional_file'])+'. Read report.md for the failed checks.')
            else:
                print('No usable correction. Read report.md for the failed check.')
            print(report.get('correction_note',''))
        return code
    except (InputError, ValueError, KeyError, OSError, EOFError) as exc:
        print(json.dumps({'status':'input_or_tool_error','error':str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 1
    except Exception as exc:
        print(json.dumps({'status':'tool_error','type':type(exc).__name__,'error':str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 1
