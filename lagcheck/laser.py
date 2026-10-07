"""Read laser intervals and export an auditable affine timestamp correction.

``header_row`` is a zero-based physical line number. CSV cells and row order are
preserved; header matching alone ignores surrounding whitespace. Numeric times
use the configured unit. Naive datetimes are calendar coordinates, never local
machine time. A sample selector is case-sensitive after trimming cell whitespace.
"""
from __future__ import annotations

import csv
import hashlib
import io
import math
import re
from datetime import datetime, timedelta, timezone
from decimal import Decimal, InvalidOperation, ROUND_HALF_EVEN
from pathlib import Path

from .models import InputError, LaserData
from . import new_wave


_UNITS = {"s": Decimal(1), "ms": Decimal("0.001"),
          "us": Decimal("0.000001"), "min": Decimal(60)}
_EPOCH = datetime(1970, 1, 1)
_MICRO = Decimal("0.000001")


def _decimal(value, name):
    try:
        result = Decimal(str(value).strip())
    except (InvalidOperation, ValueError) as exc:
        raise InputError(f"Invalid {name}: {value!r}") from exc
    if not result.is_finite():
        raise InputError(f"{name} must be finite")
    return result


def _column(header, name):
    matches = [i for i, value in enumerate(header) if value.strip() == str(name).strip()]
    if len(matches) != 1:
        raise InputError(f"Column {name!r} must occur exactly once in the header")
    return matches[0]


def _cell(row, index, row_number):
    if index >= len(row):
        raise InputError(f"Missing required cell in CSV record {row_number}")
    return row[index]


def _time(raw, config):
    value = raw.strip()
    if not value:
        raise InputError("Empty timestamp, including a metadata timestamp, is not supported")
    kind = config.get("time_kind", "datetime")
    if kind == "numeric":
        unit = config.get("time_unit")
        if unit not in _UNITS:
            raise InputError("time_unit must be s, ms, us or min")
        return _decimal(value, "timestamp") * _UNITS[unit], None
    if kind != "datetime":
        raise InputError("time_kind must be datetime or numeric")
    fmt = config.get("datetime_format")
    try:
        parsed = datetime.strptime(value, fmt) if fmt else datetime.fromisoformat(value)
    except (ValueError, TypeError) as exc:
        raise InputError(f"Cannot parse datetime {raw!r}") from exc
    # fromisoformat accepts longer fractional seconds but truncates to microseconds.
    # Retain that original fractional part exactly in the calculation.
    fraction = re.search(r"[T ]\d{2}:\d{2}:\d{2}[.,](\d+)", value) if not fmt else None
    extra = (Decimal("0." + fraction.group(1)) - Decimal(parsed.microsecond) / 1000000
             if fraction else Decimal(0))
    normalized = (parsed.astimezone(timezone.utc).replace(tzinfo=None)
                  if parsed.tzinfo is not None else parsed)
    delta = normalized - _EPOCH
    seconds = (Decimal(delta.days) * 86400 + Decimal(delta.seconds)
               + Decimal(delta.microseconds) / 1000000 + extra)
    return seconds, parsed


def _record_subset(indices, config, count):
    bounds=config.get('record_range')
    if bounds is None:
        return indices
    if (not isinstance(bounds,list) or len(bounds)!=2
            or any(isinstance(v,bool) or not isinstance(v,int) for v in bounds)
            or not 1 <= bounds[0] <= bounds[1] <= count):
        raise InputError('record_range must give two valid 1-based CSV record numbers')
    selected=[i for i in indices if bounds[0] <= i+1 <= bounds[1]]
    if not selected:
        raise InputError('record_range contains no selected records')
    return selected


def _selected_rows(rows, header, config, fmt):
    explicit = config.get("sample_column")
    sample_column = explicit if explicit is not None else ("Comment" if fmt in ("laser_state_csv", new_wave.FORMAT) else None)
    selector = config.get("sample_selector")
    if selector is not None and not isinstance(selector, str):
        raise InputError("sample_selector must be an exact sample-name string")
    if sample_column is None:
        if selector is not None:
            raise InputError("sample_selector requires sample_column")
        return _record_subset(list(range(len(rows))),config,len(rows)), [""] * len(rows), ""
    present = any(cell.strip() == str(sample_column).strip() for cell in header)
    if not present and explicit is None and selector is None:
        return _record_subset(list(range(len(rows))),config,len(rows)), [""] * len(rows), ""
    col = _column(header, sample_column)
    labels, current = [], ""
    for number, row in enumerate(rows, 1):
        raw = _cell(row, col, number).strip()
        if raw:
            current = raw
        labels.append(current)
    if config.get('sample_scope') is not None:
        if fmt != new_wave.FORMAT or config['sample_scope'] != 'all' or selector is not None:
            raise InputError('sample_scope=all is only valid for an explicit New Wave file without sample_selector')
        return _record_subset(list(range(len(rows))),config,len(rows)), labels, 'all explicitly included samples'
    names = set(labels) - {""}
    if selector is None:
        if len(names) > 1:
            raise InputError("Multiple samples: set sample_column and sample_selector")
        # With a single acquisition, retain its leading metadata as well.
        return _record_subset(list(range(len(rows))),config,len(rows)), labels, next(iter(names), "")
    indices = [i for i, label in enumerate(labels) if label == selector]
    if not indices:
        raise InputError(f"Selected sample {selector!r} was not found")
    record_range = config.get('record_range')
    if record_range is not None:
        if (not isinstance(record_range, list) or len(record_range) != 2
                or any(isinstance(v, bool) or not isinstance(v, int) for v in record_range)
                or not 1 <= record_range[0] <= record_range[1] <= len(rows)):
            raise InputError('record_range must give two valid 1-based CSV record numbers')
        indices = [i for i in indices if record_range[0] <= i + 1 <= record_range[1]]
        if not indices:
            raise InputError('record_range contains no selected sample records')
    if indices != list(range(indices[0], indices[-1] + 1)) and not config.get('_inspect_segments', False) and fmt != new_wave.FORMAT:
        raise InputError("Selected sample occurs in disjoint segments; provide one acquisition file")
    return indices, labels, selector


def _acquisition_identity(rows, header, indices, sample, config):
    """Fail before diagnosis if a selected name spans sequence identities.

    ESL blank sequence cells inherit the last explicit value. A sample label
    alone is not an acquisition identifier. Generic inputs without this field
    remain usable only under the documented pre-isolated-input contract.
    """
    present = any(cell.strip() == 'Sequence Number' for cell in header)
    expected = config.get('expected_sequence')
    if not present:
        if expected is not None:
            raise InputError('expected_sequence requires a Sequence Number column')
        return dict(sample=sample, sequence_id=None, sequence_column=None,
                    sequence_verified=False, scope='No sequence field; user must supply one isolated acquisition')
    col = _column(header, 'Sequence Number')
    labels, current = [], ''
    for number, row in enumerate(rows, 1):
        value = _cell(row, col, number).strip()
        if value:
            current = value
        labels.append(current)
    selected = [labels[i] for i in indices]
    if any(not value for value in selected):
        raise InputError('Selected acquisition has an undefined Sequence Number; provide an explicitly identified, isolated acquisition')
    sequences = sorted(set(selected))
    if len(sequences) != 1:
        raise InputError('Multiple Sequence Number values for the selected sample: '
                         + ', '.join(sequences)
                         + '. Split and verify one acquisition before diagnosis or export')
    sequence = sequences[0]
    if expected is not None and (not isinstance(expected, str) or expected != sequence):
        raise InputError('expected_sequence does not match the selected acquisition Sequence Number')
    config['expected_sequence'] = sequence
    return dict(sample=sample, sequence_id=sequence, sequence_column='Sequence Number',
                sequence_verified=True, scope='One inherited sequence identity in the selected records')


def _parse(text, config, source_path, digest, had_bom=False):
    config = dict(config)
    fmt = config.get("format", "laser_state_csv")
    if fmt not in {"laser_state_csv", "interval_csv", new_wave.FORMAT}:
        raise InputError("format must be laser_state_csv, interval_csv or new_wave_laser_log")
    delimiter = config.get("delimiter", ",")
    if not isinstance(delimiter, str) or len(delimiter) != 1:
        raise InputError("delimiter must be one character")
    header_row = config.get("header_row", 0)
    if isinstance(header_row, bool) or not isinstance(header_row, int) or header_row < 0:
        raise InputError("header_row must be a nonnegative, zero-based physical line number")
    physical = text.splitlines(keepends=True)
    if header_row >= len(physical):
        raise InputError("header_row lies beyond the input")
    preamble = "".join(physical[:header_row])
    stream = io.StringIO("".join(physical[header_row:]), newline="")
    reader = csv.reader(stream, delimiter=delimiter, strict=True)
    try:
        header = next(reader)
        header_end = reader.line_num
        rows = list(reader)
    except (csv.Error, StopIteration) as exc:
        raise InputError(f"Malformed CSV: {exc}") from exc
    if not header or not rows:
        raise InputError("Laser CSV requires a header and data records")
    header_text = "".join(physical[header_row:header_row + header_end])
    indices, labels, sample = _selected_rows(rows, header, config, fmt)
    selected = [list(rows[i]) for i in indices]
    if not selected:
        raise InputError("No selected acquisition records")
    identity = (new_wave.identity if fmt == new_wave.FORMAT else _acquisition_identity)(rows, header, indices, sample, config)
    sequences = new_wave.sequence_labels(rows, header) if fmt == new_wave.FORMAT else None
    if fmt == "interval_csv":
        if "start_column" not in config or "end_column" not in config:
            raise InputError("interval_csv requires explicit start_column and end_column")
        time_cols = [_column(header, config["start_column"]), _column(header, config["end_column"])]
        if time_cols[0] == time_cols[1]:
            raise InputError("start_column and end_column must differ")
    else:
        time_cols = [_column(header, config.get("time_column", "Timestamp"))]
    parsed_times, templates, awareness = [], [], set()
    for original, row in zip(indices, selected):
        values, row_templates = [], []
        for col in time_cols:
            try:
                value, template = _time(_cell(row, col, original + 1), config)
            except InputError as exc:
                raise InputError(f"CSV record {original + 1}: {exc}") from exc
            values.append(value)
            row_templates.append(template)
            if template is not None:
                awareness.add(template.tzinfo is not None)
        parsed_times.append(values)
        templates.append(row_templates)
    if len(awareness) > 1:
        raise InputError("Do not mix timezone-aware and naive datetimes")
    origin = parsed_times[0][0]
    intervals = []
    interval_rows = []
    if fmt == "interval_csv":
        intervals = [(values[0], values[1]) for values in parsed_times]
        interval_rows = [(i, i) for i in range(len(selected))]
    else:
        state_col = _column(header, config.get("state_column", "Laser State"))
        on = str(config.get("on_value", "On")).strip()
        off = str(config.get("off_value", "Off")).strip()
        if not on or not off or on == off:
            raise InputError("on_value and off_value must be distinct nonempty states")
        active_before = False
        for i in range(indices[0]):
            state = _cell(rows[i], state_col, i + 1).strip()
            if state == on:
                active_before = True
            elif state == off:
                active_before = False
        if active_before:
            raise InputError("Sample selection starts inside an open laser interval")
        active = None
        for i, (original, row) in enumerate(zip(indices, selected)):
            state = _cell(row, state_col, original + 1).strip()
            if active is not None and fmt == new_wave.FORMAT:
                if original != indices[i-1] + 1:
                    raise InputError('New Wave selection omits records inside a firing interval')
                if sequences[original] != sequences[indices[active]] or labels[original] != labels[indices[active]]:
                    raise InputError('New Wave firing interval crosses a line/spot sequence or sample boundary')
                if parsed_times[i][0] < parsed_times[i-1][0]:
                    raise InputError('New Wave firing state/vertex update times decrease')
            if state == on:
                if active is not None:
                    if fmt != new_wave.FORMAT:
                        raise InputError("Consecutive On records without an intervening Off")
                else:
                    active = i
            elif state == off:
                if active is not None:
                    # Prevent an unlabeled leading On from crossing a sample boundary.
                    if labels[indices[active]] != labels[original]:
                        raise InputError("Laser interval crosses a sample-name boundary")
                    intervals.append((parsed_times[active][0], parsed_times[i][0]))
                    interval_rows.append((active, i))
                    active = None
            elif state:
                raise InputError(f"Unsupported laser state {state!r} in CSV record {original + 1}")
        if active is not None:
            raise InputError("Final laser On has no Off within the selected acquisition")
    if not intervals:
        raise InputError("No complete laser firing intervals")
    policy = config.get('zero_duration_policy', 'reject')
    if policy not in ('reject', 'exclude'):
        raise InputError('zero_duration_policy must be reject or exclude')
    excluded = []
    retained, retained_rows = [], []
    previous_end = None
    for (start, end), (a, b) in zip(intervals, interval_rows):
        if end == start and policy == 'exclude':
            excluded.append(dict(reason='zero_duration', start_record=indices[a] + 1,
                                 end_record=indices[b] + 1, timestamp_seconds=str(start)))
            continue
        if end <= start:
            raise InputError("Every laser interval must have a strictly positive duration")
        if previous_end is not None and start < previous_end:
            raise InputError("Firing intervals are unordered or overlap; rows will not be sorted")
        previous_end = end
        retained.append((start, end))
        retained_rows.append((a, b))
    intervals, interval_rows = retained, retained_rows
    if not intervals:
        raise InputError('No positive-duration firing intervals remain after zero-duration exclusion')
    unordered = sum(b[0] < a[0] for a, b in zip(parsed_times, parsed_times[1:]))
    warnings = []
    if excluded:
        warnings.append(f'{len(excluded)} zero-duration firing intervals excluded from analysis; original CSV records retained unchanged')
    if unordered:
        warnings.append(f"{unordered} decreasing metadata timestamps retained in original row order; no repair applied")
    lines = [{"Start Time": float(start - origin), "End Time": float(end - origin),
              "Sample Name": sample} for start, end in intervals]
    if any(not math.isfinite(line[key]) for line in lines for key in ("Start Time", "End Time")):
        raise InputError("Relative laser times exceed the supported floating-point range")
    if any(line["End Time"] <= line["Start Time"] for line in lines):
        raise InputError("Relative firing durations are below supported floating-point precision")
    metadata = {"format": fmt, "source_path": str(source_path), "input_sha256": digest,
                "sample": sample, "selected_row_count": len(selected), "input_row_count": len(rows),
                "selected_record_start": indices[0] + 1, "selected_record_end": indices[-1] + 1,
                "time_origin": _cell(selected[0], time_cols[0], indices[0] + 1),
                "time_origin_seconds": str(origin), "time_basis": "relative seconds from first selected timestamp",
                "spot_width_um": None, "scan_speed_um_s": None, "warnings": warnings,
                "acquisition_identity": identity,
                "input_processing": dict(version='0.4.0', zero_duration_policy=policy,
                     excluded_intervals=excluded, excluded_count=len(excluded),
                     selected_record_numbers=[i + 1 for i in indices],
                     source_sha256=digest, source_path=str(source_path),
                     retained_interval_count=len(intervals),
                     note='Zero-duration pairs are excluded from analysis, never assigned a duration. Source rows are not edited.'),
                "nonmonotonic_metadata_timestamp_count": unordered,
                "config": config}
    table = {"header": header, "header_text": header_text, "preamble": preamble,
             "rows": selected, "time_columns": time_cols, "times": parsed_times,
             "templates": templates, "intervals": intervals, "interval_rows": interval_rows,
             "origin": origin, "had_bom": had_bom,
             "newline": "\r\n" if "\r\n" in text else "\n"}
    result = LaserData(lines=lines, metadata=metadata, table=table)
    return new_wave.annotate(result, indices, labels, sequences) if fmt == new_wave.FORMAT else result


def load_laser(path, config) -> LaserData:
    """Load one explicitly selected acquisition without sorting or changing cells."""
    path = Path(path)
    try:
        raw = path.read_bytes()
        text = raw.decode(config.get("encoding", "utf-8-sig"))
    except (OSError, UnicodeError, LookupError) as exc:
        raise InputError(f"Cannot read laser CSV {path}: {exc}") from exc
    return _parse(text, config, path.resolve(), hashlib.sha256(raw).hexdigest(), raw.startswith(b"\xef\xbb\xbf"))


def _output_time(value, original, template, config):
    # Six decimal places in seconds, or finer, are mandatory for rate corrections.
    if config.get("time_kind", "datetime") == "numeric":
        scaled = value / _UNITS[config.get("time_unit", "s")]
        # 12 decimal places also resolves microseconds when the unit is minutes.
        try:
            return format(scaled.quantize(Decimal("0.000000000001"), rounding=ROUND_HALF_EVEN), "f")
        except InvalidOperation as exc:
            raise InputError("Corrected numeric timestamp exceeds supported decimal precision") from exc
    micros = int((value * 1000000).to_integral_value(rounding=ROUND_HALF_EVEN))
    try:
        instant = _EPOCH + timedelta(microseconds=micros)
        if template.tzinfo is not None:
            instant = instant.replace(tzinfo=timezone.utc).astimezone(template.tzinfo)
    except (OverflowError, ValueError) as exc:
        raise InputError("Corrected datetime is outside the supported calendar range") from exc
    fmt = config.get("datetime_format")
    if fmt:
        if "%f" not in fmt:
            if "%S" not in fmt:
                raise InputError("datetime_format must contain %S to export microsecond timestamps")
            fmt = fmt.replace("%S", "%S.%f")
            config["datetime_format"] = fmt
        return instant.strftime(fmt)
    return instant.isoformat(sep="T" if "T" in original else " ", timespec="microseconds")


def export_corrected(laser: LaserData, path, rate_ppm, *, offset_seconds=0.0) -> dict:
    """Write t' = t0 + (1 + ppm/1e6) * (t - t0) + offset_seconds.

    The origin t0 is the first selected firing start. An optional offset must
    already be expressed in this file's timestamp coordinate, NOT in separately
    rebased reader coordinates. Parsing, all-cell comparison and an input hash
    check must pass before the destination is created with exclusive mode.
    """
    rate = _decimal(rate_ppm, "rate_ppm")
    offset = _decimal(offset_seconds, "offset_seconds")
    factor = Decimal(1) + rate / 1000000
    if factor <= 0:
        raise InputError("rate_ppm must produce a strictly positive time-scale factor")
    path = Path(path)
    if path.exists():
        raise InputError(f"Output already exists; refusing overwrite: {path}")
    if not isinstance(laser.table, dict) or not laser.table.get("intervals"):
        raise InputError("export_corrected requires an unmodified load_laser result")
    original_path = Path(laser.metadata["source_path"])
    original_hash = laser.metadata["input_sha256"]
    try:
        if hashlib.sha256(original_path.read_bytes()).hexdigest() != original_hash:
            raise InputError("Input file changed after loading; reload before exporting")
    except OSError as exc:
        raise InputError(f"Cannot verify original input: {exc}") from exc
    table = laser.table
    config = dict(laser.metadata["config"])
    # The output contains only selected records; original input record indices
    # cannot be reused to select rows in this shorter output.
    config.pop('record_range', None)
    config.pop('_inspect_segments', None)
    t0 = table["intervals"][0][0]
    rows = [list(row) for row in table["rows"]]
    maximum_error = Decimal(0)
    for i, row in enumerate(rows):
        for j, col in enumerate(table["time_columns"]):
            corrected = t0 + factor * (table["times"][i][j] - t0) + offset
            row[col] = _output_time(corrected, row[col], table["templates"][i][j], config)
            recovered, _ = _time(row[col], config)
            maximum_error = max(maximum_error, abs(recovered - corrected))
    buffer = io.StringIO(newline="")
    buffer.write(table["preamble"])
    buffer.write(table["header_text"])
    writer = csv.writer(buffer, delimiter=config.get("delimiter", ","), lineterminator=table["newline"])
    writer.writerows(rows)
    output_text = buffer.getvalue()
    parsed = _parse(output_text, config, path, "pending", table["had_bom"])
    if len(parsed.table["rows"]) != len(table["rows"]):
        raise InputError("Export verification failed: row count changed")
    time_cols = set(table["time_columns"])
    for before, after in zip(table["rows"], parsed.table["rows"]):
        if len(before) != len(after) or any(a != b for k, (a, b) in enumerate(zip(before, after)) if k not in time_cols):
            raise InputError("Export verification failed: non-time cells changed")
    if len(parsed.table["intervals"]) != len(table["intervals"]):
        raise InputError("Export verification failed: firing interval count changed")
    interval_error = Decimal(0)
    for before, after in zip(table["intervals"], parsed.table["intervals"]):
        for old, new in zip(before, after):
            interval_error = max(interval_error, abs(new - (t0 + factor * (old - t0) + offset)))
    if maximum_error > _MICRO or interval_error > _MICRO:
        raise InputError("Export timestamp quantization exceeded one microsecond")
    encoding = config.get("encoding", "utf-8-sig")
    if encoding.lower().replace("_", "-") == "utf-8-sig" and not table["had_bom"]:
        encoding = "utf-8"
    try:
        payload = output_text.encode(encoding)
        # The in-memory parse above verifies the exact text that will be written.
        if hashlib.sha256(original_path.read_bytes()).hexdigest() != original_hash:
            raise InputError("Input file changed during export verification")
        with path.open("xb") as destination:
            destination.write(payload)
        if path.read_bytes() != payload:
            raise InputError("Export verification failed: written bytes differ")
        unchanged = hashlib.sha256(original_path.read_bytes()).hexdigest() == original_hash
        if not unchanged:
            raise InputError("Input file changed during export")
    except (OSError, UnicodeError, LookupError) as exc:
        raise InputError(f"Cannot write corrected CSV: {exc}") from exc
    return {"output_path": str(path.resolve()), "output_sha256": hashlib.sha256(payload).hexdigest(),
            "input_sha256": original_hash, "input_unchanged": unchanged,
            "rate_ppm": float(rate), "scale_factor": str(factor),
            "origin_absolute_seconds": str(t0), "origin_relative_seconds": float(t0 - table["origin"]),
            "alignment_offset_applied": bool(offset), "offset_seconds": str(offset),
            "selected_row_count": len(rows),
            "firing_interval_count": len(table["intervals"]), "row_count_verified": True,
            "non_time_cells_verified": True, "firing_intervals_verified": True,
            "max_timestamp_quantization_error_s": float(maximum_error),
            "max_interval_endpoint_error_s": float(interval_error),
            "warnings": list(laser.metadata["warnings"]),
            "acquisition_identity": dict(laser.metadata['acquisition_identity']), "output_config": config}
