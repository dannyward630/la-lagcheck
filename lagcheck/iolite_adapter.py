"""Make an auditable millisecond, chronological ESL candidate for iolite.

This adapter neither estimates nor applies a rate/offset. Its input must already
match a completed la-lagcheck correction manifest. Actual saved-byte timing checks
and an iolite import remain separate prerequisites for accepting its output.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timedelta
from decimal import Decimal, ROUND_HALF_EVEN
import hashlib
import io
import json
from pathlib import Path
import re
import sys
from . import __version__

VERSION = __version__
ALGORITHM_VERSION = "0.2.0"
PROFILE = "esl_single_sample_chronological_milliseconds"
EPOCH = datetime(1970, 1, 1)
MILLISECOND = Decimal("0.001")
DATE_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{1,6}$")


class CompatibilityError(ValueError):
    pass


def sha256(payload):
    return hashlib.sha256(payload).hexdigest()


def timestamp_seconds(text):
    if not DATE_PATTERN.fullmatch(text):
        raise CompatibilityError("This profile requires naive ISO dates with fractional seconds and a space separator")
    try:
        instant = datetime.strptime(text, "%Y-%m-%d %H:%M:%S.%f")
    except ValueError as exc:
        raise CompatibilityError(f"Invalid timestamp: {text!r}") from exc
    delta = instant - EPOCH
    return Decimal(delta.days * 86400 + delta.seconds) + Decimal(delta.microseconds) / 1000000


def milliseconds_text(value):
    rounded = value.quantize(MILLISECOND, rounding=ROUND_HALF_EVEN)
    instant = EPOCH + timedelta(milliseconds=int(rounded * 1000))
    return instant.strftime("%Y-%m-%d %H:%M:%S.") + f"{instant.microsecond // 1000:03d}", rounded


def column(header, name):
    found = [i for i, cell in enumerate(header) if cell.strip() == name]
    if len(found) != 1:
        raise CompatibilityError(f"Required column {name!r} must occur exactly once")
    return found[0]


def inherited(rows, order, columns):
    state = [""] * len(columns)
    result = {}
    for source_id in order:
        row = rows[source_id]
        for j, col in enumerate(columns):
            if row[col].strip():
                state[j] = row[col].strip()
        result[source_id] = tuple(state)
    return result


def firing_pairs(rows, order, state_col, times, excluded_pairs=()):
    active = None
    pairs = []
    excluded_pairs = set(excluded_pairs)
    seen_excluded = set()
    last_end = None
    for source_id in order:
        state = rows[source_id][state_col].strip()
        if state == "On":
            if active is not None:
                raise CompatibilityError("Consecutive On records are unsupported")
            active = source_id
        elif state == "Off":
            if active is not None:
                start, end = times[active], times[source_id]
                if end == start and (active, source_id) in excluded_pairs:
                    seen_excluded.add((active, source_id))
                    active = None
                    continue
                if end <= start:
                    raise CompatibilityError("Firing duration is nonpositive, including after millisecond rounding")
                if last_end is not None and start < last_end:
                    raise CompatibilityError("Firing intervals overlap or are out of order")
                pairs.append((active, source_id))
                last_end = end
                active = None
        else:
            raise CompatibilityError("This ESL profile requires an explicit On or Off state in every row")
    if active is not None or not pairs:
        raise CompatibilityError("A complete nonempty set of firing intervals is required")
    if seen_excluded != excluded_pairs:
        raise CompatibilityError("Recorded zero-duration exclusions do not match actual firing pairs")
    return pairs


def prepare(payload, exclusions=()):
    try:
        text = payload.decode("utf-8-sig")
        rows = list(csv.reader(io.StringIO(text, newline=""), strict=True))
    except (UnicodeError, csv.Error) as exc:
        raise CompatibilityError(f"Cannot read UTF-8 ESL CSV: {exc}") from exc
    if len(rows) < 3:
        raise CompatibilityError("CSV needs a header and complete firing records")
    header, rows = rows[0], rows[1:]
    names = ("Timestamp", "Comment", "Sequence Number", "Scan Velocity (um/s)", "Laser State")
    tc, cc, sc, vc, lc = [column(header, name) for name in names]
    if tc != 0:
        raise CompatibilityError("This ESL profile requires Timestamp as the first column")
    if any(len(row) <= max(cc, sc, vc, lc) or len(row) > len(header) for row in rows):
        raise CompatibilityError("A CSV row lacks required fields or exceeds its header")
    if not rows[0][cc].strip() or not rows[0][sc].strip():
        raise CompatibilityError("The first record must explicitly identify its sample and sequence")
    for col, label in ((cc, "sample"), (sc, "sequence")):
        if len({row[col].strip() for row in rows if row[col].strip()}) != 1:
            raise CompatibilityError(f"Only one selected {label} is supported; split and verify acquisitions first")
    original_order = list(range(len(rows)))
    source_times = [timestamp_seconds(row[tc]) for row in rows]
    outputs = [milliseconds_text(value) for value in source_times]
    rounded_times = [value for _, value in outputs]
    excluded_pairs = [(entry['start_record'] - 1, entry['end_record'] - 1) for entry in exclusions]
    if len(set(excluded_pairs)) != len(excluded_pairs):
        raise CompatibilityError("Duplicate zero-duration exclusion")
    removed = {i for pair in excluded_pairs for i in pair}
    if any(type(i) is not int or not 0 <= i < len(rows) for i in removed):
        raise CompatibilityError("Zero-duration exclusion record is out of range")
    # Validate the complete source event stream before removing only documented zeros.
    before_pairs = firing_pairs(rows, original_order, lc, source_times, excluded_pairs)
    retained_order = [i for i in original_order if i not in removed]
    # Equal output timestamps retain their original record order.
    order = sorted(retained_order, key=lambda i: rounded_times[i])
    after_pairs = firing_pairs(rows, order, lc, rounded_times)
    if before_pairs != after_pairs:
        raise CompatibilityError("Sorting changes the identities/order of paired firing endpoints")
    before = inherited(rows, original_order, (cc, sc, vc))
    explicit_rows = [list(row) for row in rows]
    after = inherited(explicit_rows, order, (cc, sc, vc))
    materialized = []
    if any(before[i][:2] != after[i][:2] for i in order):
        for i in order:
            for j, col in enumerate((cc, sc)):
                if not explicit_rows[i][col].strip():
                    materialized.append(dict(source_record=i + 1, column=header[col],
                                             column_index=col, original=rows[i][col], value=before[i][j]))
                    explicit_rows[i][col] = before[i][j]
        after = inherited(explicit_rows, order, (cc, sc, vc))
    if any(before[i][:2] != after[i][:2] for i in order):
        raise CompatibilityError("Sorting changes inherited sample or sequence membership")
    endpoint_ids = {i for pair in before_pairs for i in pair}
    if any(before[i][2] != after[i][2] for i in endpoint_ids):
        raise CompatibilityError("Sorting changes the inherited scan velocity at a firing endpoint")
    if any(not before[i][2] for i in endpoint_ids):
        raise CompatibilityError("Every firing endpoint requires an explicit preceding scan velocity")
    # Only recorded empty membership cells can be materialized; all other cells
    # (including coordinates and scan velocity) remain byte-for-cell equal.
    output_rows = []
    mapping = []
    starts = {pair[0] for pair in before_pairs}
    ends = {pair[1] for pair in before_pairs}
    for position, i in enumerate(order, 1):
        row = list(explicit_rows[i])
        row[tc] = outputs[i][0]
        output_rows.append(row)
        mapping.append(dict(output_record=position, source_record=i + 1,
                            source_timestamp=rows[i][tc], output_timestamp=row[tc],
                            quantization_error_s=str(rounded_times[i] - source_times[i]),
                            role="firing_start" if i in starts else "firing_end" if i in ends else "metadata",
                            sample=before[i][0], sequence=before[i][1],
                            inherited_speed_before=before[i][2], inherited_speed_after=after[i][2]))
    stream = io.StringIO(newline="")
    newline = "\r\n" if "\r\n" in text else "\n"
    writer = csv.writer(stream, lineterminator=newline)
    writer.writerow(header)
    writer.writerows(output_rows)
    output = stream.getvalue().encode("utf-8-sig" if payload.startswith(b"\xef\xbb\xbf") else "utf-8")
    recovered = list(csv.reader(io.StringIO(output.decode("utf-8-sig"), newline=""), strict=True))
    if recovered[0] != header or recovered[1:] != output_rows:
        raise CompatibilityError("Output CSV did not round-trip exactly")
    if any(recovered[position][1:] != explicit_rows[i][1:] for position, i in enumerate(order, 1)):
        raise CompatibilityError("A mapped non-time cell changed")
    errors = [abs(rounded_times[i] - source_times[i]) for i in original_order]
    if max(errors) > Decimal("0.0005"):
        raise CompatibilityError("Millisecond rounding exceeded half a millisecond")
    details = dict(profile=PROFILE, selected_sample=before[0][0], selected_sequence=before[0][1],
                   input_row_count=len(rows), output_row_count=len(output_rows),
                   firing_interval_count=len(before_pairs), moved_row_count=sum(i != p for p, i in enumerate(order)),
                   timestamp_precision="milliseconds, exactly three fractional digits",
                   rounding="ROUND_HALF_EVEN", sorting="stable by rounded timestamp, original order for ties",
                   timestamp_interpolation_applied=False, rows_removed=len(removed), rate_or_offset_applied=False,
                   excluded_zero_duration_intervals=list(exclusions),
                   excluded_source_records=[dict(source_record=i + 1, cells=rows[i]) for i in sorted(removed)],
                   materialized_membership_cells=materialized,
                   max_timestamp_quantization_error_s=str(max(errors)),
                   non_time_cells_equal_by_source_record=not materialized,
                   non_membership_cells_equal_by_source_record=True,
                   sample_and_sequence_inheritance_equal_all_records=True,
                   scan_velocity_inheritance_equal_firing_endpoints=True,
                   metadata_records_with_changed_inherited_speed=[i + 1 for i in order
                       if i not in endpoint_ids and before[i][2] != after[i][2]],
                   firing_endpoint_source_pairs_equal=True,
                   firing_endpoint_source_records=[[a + 1, b + 1] for a, b in before_pairs],
                   source_record_by_output_record=[i + 1 for i in order])
    return output, mapping, details


def write_json(path, data):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def recorded_exclusions(canonical, provenance, payload):
    """Bind exclusion permission to the completed report and unchanged raw input.

    A zero in the corrected CSV alone is never permission to remove an event.
    Original record numbers are translated through the recorded selection map.
    """
    report_path = canonical.parent / 'report.json'
    if not report_path.exists():
        return [], {}
    report_bytes = report_path.read_bytes()
    report = json.loads(report_bytes.decode('utf-8-sig'))
    metadata = (report.get('input_metadata') or {}).get('laser') or {}
    processing = metadata.get('input_processing') or {}
    excluded = processing.get('excluded_intervals') or []
    if not excluded:
        return [], {}
    config_path = canonical.parent / 'config.json'
    config_bytes = config_path.read_bytes()
    config = json.loads(config_bytes.decode('utf-8-sig'))
    options = config.get('laser') or {}
    if (processing.get('zero_duration_policy') != 'exclude'
            or options.get('zero_duration_policy') != 'exclude'
            or processing.get('excluded_count') != len(excluded)):
        raise CompatibilityError('Zero-duration exclusion lacks an explicit matching policy/count')
    receipt_key = 'provisional_export_verification' if provenance.get('provisional') else 'export_verification'
    receipt = report.get(receipt_key) or {}
    expected_hash = (provenance.get('input_sha256') or {}).get('laser')
    if (not expected_hash or report.get('input_sha256') != provenance.get('input_sha256')
            or metadata.get('input_sha256') != expected_hash
            or processing.get('source_sha256') != expected_hash
            or receipt.get('input_sha256') != expected_hash
            or receipt.get('output_sha256') != sha256(payload)):
        raise CompatibilityError('Zero-duration exclusion provenance hashes do not match')
    raw_path = (canonical.parent / options['path']).resolve()
    if sha256(raw_path.read_bytes()) != expected_hash:
        raise CompatibilityError('Original laser input changed since zero-duration exclusion')
    from .laser import load_laser
    laser = load_laser(raw_path, {k: v for k, v in options.items() if k != 'path'})
    fresh = laser.metadata.get('input_processing') or {}
    for key in ('zero_duration_policy', 'excluded_intervals', 'excluded_count',
                'selected_record_numbers', 'source_sha256', 'retained_interval_count'):
        if fresh.get(key) != processing.get(key):
            raise CompatibilityError('Recorded zero-duration exclusion disagrees with original input: ' + key)
    source_rows = list(csv.reader(io.StringIO(payload.decode('utf-8-sig'), newline=''), strict=True))[1:]
    records = processing['selected_record_numbers']
    if (len(records) != len(source_rows) or len(set(records)) != len(records)
            or len(laser.table['rows']) != len(source_rows)):
        raise CompatibilityError('Exclusion selection map does not match the canonical CSV')
    scale = Decimal(str(receipt['scale_factor']))
    anchor = Decimal(str(receipt['origin_absolute_seconds']))
    offset = Decimal(str(receipt['offset_seconds']))
    if (scale != Decimal(1) + Decimal(str(provenance['rate_ppm'])) / 1000000
            or offset != Decimal(str(provenance['export_plan']['offset_seconds']))):
        raise CompatibilityError('Exclusion receipt disagrees with the canonical time mapping')
    for old, new, times in zip(laser.table['rows'], source_rows, laser.table['times']):
        if new[1:] != old[1:]:
            raise CompatibilityError('Exclusion map changes original non-time cells')
        expected = anchor + scale * (times[0] - anchor) + offset
        if abs(timestamp_seconds(new[0]) - expected) > Decimal('0.000001'):
            raise CompatibilityError('Exclusion map disagrees with corrected source timestamps')
    local = {record: i + 1 for i, record in enumerate(records)}
    result = []
    for entry in excluded:
        if entry.get('reason') != 'zero_duration':
            raise CompatibilityError('Only recorded zero-duration exclusions are supported')
        a, b = entry['start_record'], entry['end_record']
        if a not in local or b not in local:
            raise CompatibilityError('Excluded original record is outside the selected acquisition')
        result.append(dict(reason='zero_duration', start_record=local[a], end_record=local[b],
                           original_start_record=a, original_end_record=b,
                           original_timestamp_seconds=entry['timestamp_seconds']))
    snapshots = {str(report_path): sha256(report_bytes), str(config_path): sha256(config_bytes),
                 str(raw_path): expected_hash}
    return result, snapshots


def adapt(canonical, manifest, output_directory, allow_provisional=False):
    canonical, manifest, output_directory = map(lambda p: Path(p).resolve(), (canonical, manifest, output_directory))
    if output_directory.exists():
        raise CompatibilityError("Output directory already exists; no files are overwritten")
    if (canonical.parent / "INCOMPLETE.json").exists() or (manifest.parent / "INCOMPLETE.json").exists():
        raise CompatibilityError("An INCOMPLETE source run cannot be adapted")
    payload, manifest_bytes = canonical.read_bytes(), manifest.read_bytes()
    source_hash, manifest_hash = sha256(payload), sha256(manifest_bytes)
    try:
        provenance = json.loads(manifest_bytes.decode("utf-8-sig"))
    except (ValueError, UnicodeError) as exc:
        raise CompatibilityError("Cannot parse source correction manifest") from exc
    if not isinstance(provenance, dict) or provenance.get("output_sha256") != source_hash:
        raise CompatibilityError("Canonical CSV does not match the source manifest hash")
    provisional_source = allow_provisional and provenance.get("provisional") is True and provenance.get("verified") is False
    if provenance.get("compatibility_adapter") or (provenance.get("full_export_diagnostic_passed") is not True
                                                   and not provisional_source):
        raise CompatibilityError("Input must be a completed canonical correction, not an already adapted file")
    if provenance.get("mode") not in ("affine", "rate_only"):
        raise CompatibilityError("Source manifest must identify affine or rate_only mode")
    exclusions, exclusion_sources = recorded_exclusions(canonical, provenance, payload)
    transformed, mapping, details = prepare(payload, exclusions) if exclusions else prepare(payload)
    details['exclusion_source_sha256'] = exclusion_sources
    def check_exclusion_sources():
        if any(sha256(Path(path).read_bytes()) != expected for path, expected in exclusion_sources.items()):
            raise CompatibilityError('Zero-duration exclusion evidence changed during adaptation')
    check_exclusion_sources()
    if canonical.read_bytes() != payload or manifest.read_bytes() != manifest_bytes:
        raise CompatibilityError("Input CSV or manifest changed during preparation")
    output_directory.mkdir(parents=True, exist_ok=False)
    try:
        target = output_directory / "corrected_laser_log_iolite.csv"
        with target.open("xb") as stream:
            stream.write(transformed)
        saved = target.read_bytes()
        if saved != transformed:
            raise CompatibilityError("Saved compatibility CSV bytes differ from the prepared bytes")
        # Re-parse the actual saved bytes, then verify all mapped cells/dates.
        saved_rows = list(csv.reader(io.StringIO(saved.decode("utf-8-sig"), newline=""), strict=True))[1:]
        source_rows = list(csv.reader(io.StringIO(payload.decode("utf-8-sig"), newline=""), strict=True))[1:]
        materialized = {(entry['source_record'], entry['column_index']): entry
                        for entry in details['materialized_membership_cells']}
        for row, entry in zip(saved_rows, mapping):
            expected_row = list(source_rows[entry['source_record'] - 1])
            for col in range(1, len(expected_row)):
                change = materialized.get((entry['source_record'], col))
                if change:
                    if expected_row[col] != change['original'] or expected_row[col].strip():
                        raise CompatibilityError('Only empty inherited membership cells may be materialized')
                    expected_row[col] = change['value']
            if row[1:] != expected_row[1:] or row[0] != entry["output_timestamp"]:
                raise CompatibilityError("Saved compatibility CSV fails its source-record mapping")
        if len(saved_rows) != len(source_rows) - details['rows_removed'] or len(saved_rows) != len(mapping):
            raise CompatibilityError("Saved compatibility CSV changed the row count")
        if details['excluded_source_records']:
            write_json(output_directory / 'excluded_zero_duration_records.json',
                       dict(intervals=exclusions, records=details['excluded_source_records'],
                            canonical_source_sha256=source_hash, evidence_sha256=exclusion_sources,
                            note='Only already excluded zero-duration On/Off records are omitted from the compatibility CSV; their exact cells and source identities are retained here.'))
        with (output_directory / "row_mapping.csv").open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(mapping[0]))
            writer.writeheader()
            writer.writerows(mapping)
        output_hash = sha256(saved)
        receipt = dict(adapter_version=VERSION, adapter_algorithm_version=ALGORITHM_VERSION, canonical_source_path=str(canonical),
                       canonical_source_sha256=source_hash, canonical_manifest_path=str(manifest),
                       canonical_manifest_sha256=manifest_hash, canonical_mode=provenance["mode"],
                       output_path=str(target), output_sha256=output_hash, **details,
                       saved_bytes_match_prepared=True, full_export_diagnostic_passed=False,
                       importer_roundtrip_tested=False, usable_file=None,
                       status="compatibility_candidate_requires_validation",
                       required_followup="Run the full timing diagnostic on these saved bytes with the declared fixed mapping, then verify this exact hash in the receiving iolite version. Do not apply a rate or original offset again.")
        if canonical.read_bytes() != payload or manifest.read_bytes() != manifest_bytes:
            raise CompatibilityError("Input CSV or manifest changed during output")
        check_exclusion_sources()
        write_json(output_directory / "compatibility_receipt.json", receipt)
        # This is the last marker. The legacy CLI recognizes output_sha256 and
        # consequently refuses to apply another correction to this same file.
        write_json(output_directory / "correction_manifest.json",
                   dict(version=f"iolite-compatibility-adapter-{VERSION}", compatibility_adapter=True,
                        profile=PROFILE, mode=provenance["mode"], output_sha256=output_hash,
                        output_file=str(target), canonical_source_sha256=source_hash,
                        canonical_manifest_sha256=manifest_hash,
                        canonical_export_plan=provenance.get("export_plan"),
                        canonical_tool_version=provenance.get("version"),
                        full_export_diagnostic_passed=False, importer_roundtrip_tested=False,
                        status=receipt["status"], instruction=receipt["required_followup"]))
        return receipt
    except Exception as exc:
        write_json(output_directory / "INCOMPLETE.json", dict(error=str(exc), usable_file=None,
                   instruction="Do not use any files from this incomplete adaptation"))
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("canonical_csv")
    parser.add_argument("canonical_manifest")
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        receipt = adapt(args.canonical_csv, args.canonical_manifest, args.output)
    except (CompatibilityError, OSError, ValueError) as exc:
        print(json.dumps(dict(status="rejected", error=str(exc))), file=sys.stderr)
        return 1
    print(json.dumps({k: receipt[k] for k in ("status", "output_path", "output_sha256",
                     "firing_interval_count", "max_timestamp_quantization_error_s")}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
