"""Two-stage iolite export with the canonical run's fixed timing comparison."""
from __future__ import annotations

import json
from pathlib import Path

from . import __version__
from .iolite_adapter import adapt
from .laser import load_laser
from .models import InputError
from .ms import load_signal
from .verification import verify_written_log


def export_iolite(canonical_run, output):
    """Adapt a completed run, diagnose saved bytes, then publish a new file.

    This does not invoke iolite, reapply a correction, or grant importer-tested
    status. Existing canonical reports are required so no time origin or
    comparison offset is inferred from the adapted file. A run that wrote only a
    provisional (not verified) log is formatted the same way; its output stays
    marked as provisional.
    """
    from .cli import digest, write_json, write_blocks

    origin, out = Path(canonical_run).resolve(), Path(output).resolve()
    if out.exists():
        raise InputError('Output directory must not already exist; no files are overwritten')
    if (origin / 'INCOMPLETE.json').exists():
        raise InputError('An INCOMPLETE canonical run cannot be exported')
    provisional = (not (origin / 'correction_manifest.json').exists()
                   and (origin / 'provisional_manifest.json').exists())
    csv_name = 'provisional_laser_log.csv' if provisional else 'corrected_laser_log.csv'
    manifest_name = 'provisional_manifest.json' if provisional else 'correction_manifest.json'
    files = {name: origin / name for name in ('config.json', 'report.json', manifest_name, csv_name)}
    snapshots = {str(path): digest(path) for path in files.values()}
    config, report, manifest = [json.loads(files[name].read_text(encoding='utf-8-sig'))
                                for name in ('config.json', 'report.json', manifest_name)]
    if any(not isinstance(item, dict) for item in (config, report, manifest)):
        raise InputError('Canonical config, report and manifest must be objects')
    if manifest.get('compatibility_adapter'):
        raise InputError('Use the canonical run, not an already adapted output')
    if provisional:
        if (manifest.get('provisional') is not True or manifest.get('verified') is not False
                or report.get('correction_exported') is not False
                or report.get('provisional_file') is None):
            raise InputError('Provisional run records are inconsistent')
        report_plan = report.get('provisional_export_plan')
        source_receipt = report.get('provisional_export_verification', {})
    else:
        if (manifest.get('full_export_diagnostic_passed') is not True
                or report.get('correction_exported') is not True
                or report.get('status') != 'supported'
                or report.get('post_export_verification', {}).get('passed') is not True):
            raise InputError('Canonical correction and actual-file diagnostic must both have passed')
        report_plan = report.get('export_plan')
        source_receipt = report.get('export_verification', {})
    if manifest.get('output_sha256') != snapshots[str(files[csv_name])]:
        raise InputError('Canonical CSV does not match the completed manifest')
    plan = manifest.get('export_plan')
    if (not isinstance(plan, dict) or plan != report_plan
            or plan.get('mode') != manifest.get('mode')
            or plan.get('mode') not in ('rate_only', 'affine')):
        raise InputError('Canonical report and manifest disagree about the fixed export mapping')
    if (source_receipt.get('output_sha256') != manifest['output_sha256']
            or not isinstance(source_receipt.get('output_config'), dict)):
        raise InputError('Canonical export receipt does not describe the completed CSV')
    expected_inputs = manifest.get('input_sha256')
    if not isinstance(expected_inputs, dict) or expected_inputs != report.get('input_sha256'):
        raise InputError('Canonical report and manifest disagree about input hashes')
    paths = {}
    for kind in ('ms', 'laser'):
        section = config.get(kind)
        if not isinstance(section, dict) or not isinstance(section.get('path'), str):
            raise InputError('Canonical config requires explicit ms and laser paths')
        paths[kind] = (origin / section['path']).resolve()
        actual = digest(paths[kind])
        if actual != expected_inputs.get(kind):
            raise InputError(kind + ' input changed since the canonical run')
        snapshots[str(paths[kind])] = actual

    def check_sources():
        if any(digest(path) != expected for path, expected in snapshots.items()):
            raise InputError('A canonical record or original input changed during iolite export')

    ms_options = {k: v for k, v in config['ms'].items() if k != 'path'}
    if ms_options.get('format') == 'netcdf':
        from .netcdf import load_netcdf_signal
        ms = load_netcdf_signal(paths['ms'], ms_options)
    else:
        ms = load_signal(paths['ms'], ms_options)
    laser = load_laser(paths['laser'], {k: v for k, v in config['laser'].items() if k != 'path'})
    check_sources()
    out.mkdir(parents=True, exist_ok=False)
    try:
        compatibility = adapt(files[csv_name], files[manifest_name], out / 'format_stage',
                              allow_provisional=provisional)
        candidate = Path(compatibility['output_path'])
        receipt = dict(source_receipt, output_path=str(candidate), output_sha256=digest(candidate))
        # The canonical output reader accepts fractional ISO dates. Sorting
        # changes record positions, so the adapter separately validates every
        # endpoint identity and each non-time cell by its source record.
        post = None
        fallback = report.get('alignment_fallback') or {}
        timing_report = fallback['diagnostic'] if provisional and fallback.get('adopted') else report
        if (timing_report.get('alignment') or {}).get('offset_s') is not None:
            post, rerun, blocks = verify_written_log(
                ms, laser, timing_report, candidate, receipt, plan, config['analysis'])
            write_json(out / 'post_export_report.json', rerun)
            write_blocks(out, blocks, 'post_export_blocks')
            write_json(out / 'iolite_timing_verification.json', post)
        elif not provisional:
            raise InputError('Canonical report lacks the alignment needed for verification')
        instruction = ('Affine timestamps already contain the fitted rate and constant. '
                       'For a matching verified iolite coordinate, disable Auto sync and use Laser Offset 0 s; '
                       'verify this file in the actual importer before accepting the workspace.'
                       if plan['mode'] == 'affine' else
                       'The rate is already corrected. Determine only the constant offset in iolite. '
                       'Before image use, verify that firing-interval endpoints and continuous-trajectory transition times '
                       'use the same corrected time coordinate after the offset. '
                       'Do not copy another acquisition offset or reapply the rate.')
        if provisional:
            instruction = ('PROVISIONAL, NOT VERIFIED: the timing checks of the canonical run did not all pass (stage: '
                           + str(manifest.get('failure_stage') or 'input_diagnostic') + '; '
                           + (', '.join(manifest.get('failed_checks') or []) or str(manifest.get('reason') or 'see report'))
                           + '). Inspect the reconstruction and state the failed checks with any result. ' + instruction)
            if fallback.get('adopted'):
                instruction += ' ' + fallback['note']
        result = dict(tool='la-lagcheck', version=__version__, canonical_run=str(origin),
                      canonical_mode=plan['mode'], canonical_input_sha256=expected_inputs,
                      canonical_records_sha256={Path(p).name: h for p, h in snapshots.items()
                                                if Path(p).parent == origin},
                      compatibility=compatibility, post_export_verification=post,
                      rate_or_offset_applied=False, correction_exported=False, usable_file=None,
                      provisional=provisional, provisional_file=None,
                      alignment_fallback=fallback,
                      importer_roundtrip_tested=False, gui_workflow_tested=False,
                      acquisition_identity=laser.metadata['acquisition_identity'],
                      receiver_interval_trajectory_agreement_verified=False,
                      independent_geometry_tested=False,
                      status='compatibility_export_qc_failed',
                      receiving_software_instruction=instruction,
                      scope='Internal checks of this exact output, not a new iolite application test')
        check_sources()
        verified = post is not None and post['passed'] and not provisional
        if verified or provisional:
            # Bind publication to the bytes the adapter checked, whether or not a timing read-back ran.
            expected = compatibility['output_sha256']
            payload = candidate.read_bytes()
            if digest(candidate) != expected:
                raise InputError('Compatibility candidate changed relative to the checked adapter output')
            if post is not None and expected != post['diagnostic_input_sha256']:
                raise InputError('Compatibility candidate changed after the timing diagnostic')
            final = out / ('provisional_laser_log_iolite.csv' if provisional else 'corrected_laser_log_iolite.csv')
            with final.open('xb') as stream:
                stream.write(payload)
            if digest(final) != expected:
                raise InputError('Published CSV differs from the diagnosed compatibility bytes')
            if provisional:
                result.update(status='provisional_not_verified', provisional_file=str(final),
                              output_sha256=digest(final))
            else:
                result.update(status='internal_checks_passed_importer_validation_pending',
                              correction_exported=True, usable_file=str(final), output_sha256=digest(final))
        check_sources()
        write_json(out / 'iolite_export_report.json', result)
        (out / 'report.md').write_text('\n'.join([
            '# iolite compatibility export', '', 'Status: ' + result['status'],
            'Version: ' + __version__, 'Usable file: ' + str(result['usable_file']),
            'Provisional (not verified) file: ' + str(result['provisional_file']),
            'Full saved-byte timing checks: ' + str(None if post is None else post['passed']),
            'Actual importer / GUI tested by this command: False', '',
            result['receiving_software_instruction'], '',
            'format_stage contains the preserved formatting-only candidate and its row map.',
            'Only a final correction_manifest.json without INCOMPLETE.json marks a verified publication;',
            'a provisional_manifest.json marks a provisional file.',
            'Existing project validation applies only to its recorded hashes and settings.', '']), encoding='utf-8')
        if result['correction_exported'] or result['provisional_file']:
            check_sources()
            name = 'provisional_manifest.json' if provisional else 'correction_manifest.json'
            carried = ({k: manifest.get(k) for k in ('rate_ppm', 'rate_source', 'diagnostic_status', 'failed_checks',
                                                    'failure_stage', 'readback_status', 'alignment_fallback')} if provisional else {})
            if provisional:
                carried['iolite_readback_status'] = ('not_run' if post is None else
                                                     'passed' if post['passed'] else 'timing_checks_failed')
            write_json(out / name, dict(
                carried, version=__version__, compatibility_adapter=True, mode=plan['mode'],
                provisional=provisional, verified=not provisional,
                output_file=result['usable_file'] or result['provisional_file'], output_sha256=result['output_sha256'],
                canonical_output_sha256=manifest['output_sha256'],
                canonical_export_plan=plan, input_sha256=expected_inputs,
                full_export_diagnostic_passed=not provisional, importer_roundtrip_tested=False,
                acquisition_identity=laser.metadata['acquisition_identity'],
                receiver_interval_trajectory_agreement_verified=False,
                gui_workflow_tested=False, status=result['status'],
                instruction=result['receiving_software_instruction']))
        return (0 if verified else 2), result
    except Exception as exc:
        write_json(out / 'INCOMPLETE.json', dict(status='tool_or_input_error', error=str(exc),
                   usable_file=None, instruction='Do not use outputs from this incomplete run'))
        raise
