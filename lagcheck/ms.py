"""Strict signal import with explicit channel selection and time provenance.

The named MassHunter adapter recognizes its ``Intensity Vs Time`` preamble,
``Time [Sec]`` header, and optional ``Printed:`` footer. The Qtegra adapter
recognizes the iCAP RQ LDR CSV layout, including its channel metadata row.
Generic CSV imports require a column mapping; they do not imply validation
of any other vendor.
"""

import csv
from datetime import datetime
import hashlib
import io
from pathlib import Path
import re
from collections.abc import Mapping

import numpy as np

from .models import InputError, SignalData


_UNIT_SECONDS = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "min": 60.0}
_ISOTOPE = re.compile(r"^[A-Z][a-z]?\d+(?:\.[0-9]+)?$")
_QTEGRA_ISOTOPE = re.compile(r"^\d+[A-Z][a-z]?$")


def _qtegra_preamble(rows):
    """Recognize only the observed File Version 1, iCAP RQ LDR layout."""
    if len(rows) != 13 or any(len(row) != 1 for row in rows[:12]) or rows[12]:
        raise InputError("Qtegra LDR requires its 12 metadata records and blank separator.")
    identifier = re.fullmatch(
        r"([^:\r\n]+):(\d{2}-\d{2}-\d{4} \d{2}:\d{2}:\d{2} [AP]M);", rows[0][0])
    software = re.fullmatch(
        r"Software:Name=Qtegra;Version=(\d+(?:\.\d+)+);File Version=1;", rows[1][0])
    if identifier is None or software is None:
        raise InputError("Qtegra LDR identifier/Software/version metadata is missing or unsupported.")
    try:
        datetime.strptime(identifier.group(2), "%m-%d-%Y %I:%M:%S %p")
    except ValueError as exc:
        raise InputError("Qtegra LDR identifier timestamp is invalid.") from exc
    if rows[2][0] != "Configuration:Machine=iCAP RQ;":
        raise InputError("This Qtegra LDR adapter recognizes only the observed iCAP RQ layout.")
    prefixes = ("STD:", "RF Generator:", "Ion Optics:", "Vacuum:", "Detector:",
                "Cooling System:", "Power Supply:", "Gas Supply:", "Pulse Counting:")
    if any(not row[0].startswith(prefix) or not row[0].endswith(";")
           for row, prefix in zip(rows[3:12], prefixes)):
        raise InputError("Qtegra LDR instrument metadata layout is not recognized.")
    return {"software": "Qtegra", "software_version": software.group(1),
            "file_version": "1", "machine": "iCAP RQ",
            "sample_identifier": identifier.group(1),
            "identifier_timestamp_raw": identifier.group(2)}


def _qtegra_channel_metadata(row, header, row_number):
    if len(row) != len(header) + 1 or row[0].strip() or row[-1].strip():
        raise InputError(f"Row {row_number}: invalid Qtegra LDR channel metadata row.")
    metadata = {}
    for name, field in zip(header[1:], row[1:-1]):
        match = re.fullmatch(r"dwell time=([^;]+);xcal factor=([^;]+)", field.strip())
        if match is None:
            raise InputError(f"Row {row_number}: missing Qtegra dwell time/xcal metadata for {name}.")
        dwell = _number(match.group(1), f"{name} dwell time", row_number)
        xcal = _number(match.group(2), f"{name} xcal factor", row_number)
        if dwell <= 0 or xcal <= 0:
            raise InputError(f"Row {row_number}: dwell time and xcal factor must be positive.")
        # The export does not label these units. Preserve the field without inferring units.
        metadata[name] = {"dwell_time_raw": match.group(1), "xcal_factor": xcal}
    return metadata


def _required(config, name):
    if name not in config:
        raise InputError(f"Signal configuration requires {name!r}.")
    return config[name]


def _channels(config):
    channels = _required(config, "signal_columns")
    if (not isinstance(channels, list) or not channels
            or any(not isinstance(c, str) or not c.strip() for c in channels)):
        raise InputError("signal_columns must be a nonempty list of column names.")
    if len(set(channels)) != len(channels):
        raise InputError("signal_columns contains duplicate names.")
    return channels


def _number(value, column, row_number):
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise InputError(
            f"Row {row_number}: column {column!r} is not a finite number: {value!r}."
        ) from exc
    if not np.isfinite(number):
        raise InputError(f"Row {row_number}: column {column!r} is not finite.")
    return number


def _datetime(value, date_format, row_number):
    value = value.strip()
    try:
        if date_format is not None:
            return datetime.strptime(value, date_format)
        if not re.match(r"^\d{4}-\d{2}-\d{2}[Tt ]\d{2}:\d{2}", value):
            raise ValueError("ISO 8601 date and time required")
        return datetime.fromisoformat(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise InputError(
            f"Row {row_number}: invalid datetime {value!r}; provide "
            "datetime_format for a non-ISO date/time."
        ) from exc


def load_signal(path, config):
    """Return selected-channel sum with a strictly increasing time axis in seconds.

    All formats require ``signal_columns``. For ``agilent_masshunter_csv``,
    the exact named header supplies the numeric seconds mapping; ``encoding``
    defaults to UTF-8 with an optional BOM. ``thermo_qtegra_ldr_csv`` recognizes
    the File Version 1 iCAP RQ LDR layout with one Time column, mass-first
    isotope names, a dwell/xcal row and a trailing empty field. It requires
    explicit ``time_unit`` because the header does not state a unit. Its
    observed-format fixture has only five points, Qtegra 2.7.2425.65; this is
    not long-run validation or support for ordinary per-channel Time exports.
    For ``csv``, explicitly provide
    ``header_row`` (zero-based physical line), ``delimiter``, ``encoding``,
    ``time_column``, and ``time_kind`` (``numeric`` or ``datetime``). Numeric
    time also requires ``time_unit`` in s, ms, us, or min. Datetime accepts ISO
    8601 or an explicit Python ``datetime_format`` and requires consistent
    timezone awareness. Generic preamble text is sliced off before CSV parsing;
    quoted line breaks in it do not change the configured physical line number.
    No missing data are dropped, sorted, or interpolated.
    """
    if not isinstance(config, Mapping):
        raise InputError("Signal configuration must be a mapping.")
    source_format = _required(config, "format")
    if source_format not in ("agilent_masshunter_csv", "thermo_qtegra_ldr_csv", "csv"):
        raise InputError(f"Unsupported signal format: {source_format!r}.")
    channels = _channels(config)
    named = source_format == "agilent_masshunter_csv"
    qtegra = source_format == "thermo_qtegra_ldr_csv"

    if named:
        delimiter = config.get("delimiter", ",")
        if delimiter != ",":
            raise InputError("agilent_masshunter_csv requires a comma delimiter.")
        encoding = config.get("encoding", "utf-8-sig")
        for key, expected in (("time_kind", "numeric"), ("time_unit", "s"),
                              ("time_column", "Time [Sec]")):
            if key in config and config[key] != expected:
                raise InputError(f"agilent_masshunter_csv requires {key}={expected!r}.")
        time_kind, time_unit, time_column = "numeric", "s", "Time [Sec]"
        header_row = None
    elif qtegra:
        delimiter = config.get("delimiter", ",")
        if delimiter != ",":
            raise InputError("thermo_qtegra_ldr_csv requires a comma delimiter.")
        encoding = config.get("encoding", "utf-8-sig")
        for key, expected in (("time_kind", "numeric"), ("time_column", "Time")):
            if key in config and config[key] != expected:
                raise InputError(f"thermo_qtegra_ldr_csv requires {key}={expected!r}.")
        if "header_row" in config and (type(config["header_row"]) is not int
                                       or config["header_row"] != 13):
            raise InputError("thermo_qtegra_ldr_csv requires header_row=13.")
        time_kind, time_column, header_row = "numeric", "Time", 13
        time_unit = _required(config, "time_unit")
    else:
        header_row = _required(config, "header_row")
        if (isinstance(header_row, bool) or not isinstance(header_row, int)
                or header_row < 0):
            raise InputError("header_row must be a nonnegative, zero-based physical line number.")
        delimiter = _required(config, "delimiter")
        encoding = _required(config, "encoding")
        time_column = _required(config, "time_column")
        time_kind = _required(config, "time_kind")
        if not isinstance(time_column, str) or not time_column:
            raise InputError("time_column must be a column name.")
        if time_kind not in ("numeric", "datetime"):
            raise InputError("time_kind must be 'numeric' or 'datetime'.")
        time_unit = _required(config, "time_unit") if time_kind == "numeric" else None
    if (not isinstance(delimiter, str) or len(delimiter) != 1
            or delimiter in "\r\n\x00"):
        raise InputError("delimiter must be one non-newline character.")
    if not isinstance(encoding, str) or not encoding:
        raise InputError("encoding must be an explicit codec name.")
    if time_kind == "numeric" and (
            not isinstance(time_unit, str) or time_unit not in _UNIT_SECONDS):
        raise InputError(f"Unsupported time_unit {time_unit!r}; use s, ms, us, or min.")
    if time_column in channels:
        raise InputError("The time column cannot also be a signal column.")
    date_format = config.get("datetime_format")
    if date_format is not None and (not isinstance(date_format, str) or not date_format):
        raise InputError("datetime_format must be a nonempty strptime format string.")

    try:
        source = Path(path).resolve()
        raw = source.read_bytes()
        text = raw.decode(encoding)
    except (OSError, UnicodeError, LookupError, TypeError, ValueError) as exc:
        raise InputError(f"Cannot read signal file {path!s}: {exc}") from exc
    digest = hashlib.sha256(raw).hexdigest()
    physical_offset = 0
    preamble_physical_text = None
    if source_format == "csv":
        physical = text.splitlines(keepends=True)
        if header_row >= len(physical):
            raise InputError("header_row lies beyond the signal input.")
        physical_offset = header_row
        preamble_physical_text = "".join(physical[:header_row])
        text = "".join(physical[header_row:])
    reader = csv.reader(io.StringIO(text, newline=""), delimiter=delimiter, strict=True)
    header = None
    intensity_marker = False
    preamble = []
    footer = []
    time_values, signal_values = [], []
    first_value = None
    first_raw = None
    pending_blank = 0
    trailer_blanks = 0
    in_footer = False
    qtegra_metadata = {}
    awaiting_qtegra_metadata = False
    try:
        pending_delimited = []
        footer_padding_records = []
        for record_index, row in enumerate(reader):
            row_number = physical_offset + reader.line_num
            if header is None:
                if named:
                    if row and row[0].strip() == "Intensity Vs Time":
                        intensity_marker = True
                    found = bool(row) and row[0].strip() == "Time [Sec]"
                elif qtegra:
                    found = record_index == header_row
                else:
                    found = record_index == 0
                if not found:
                    preamble.append(row)
                    continue
                if named and not intensity_marker:
                    raise InputError(
                        "Time [Sec] header found without the MassHunter "
                        "Intensity Vs Time preamble; use an explicit generic CSV mapping."
                    )
                if qtegra:
                    qtegra_metadata = _qtegra_preamble(preamble)
                    if (len(row) < 3 or row[0].strip() != "Time"
                            or row[-1].strip()):
                        raise InputError("Qtegra LDR requires a Time-first header and trailing empty field.")
                    row = row[:-1]
                header = [value.strip() for value in row]
                header_row = physical_offset if source_format == "csv" else record_index
                if any(not name for name in header) or len(set(header)) != len(header):
                    raise InputError("CSV header contains empty or duplicate column names.")
                missing = [name for name in [time_column] + channels if name not in header]
                if missing:
                    raise InputError(f"Columns not found in signal file: {missing!r}.")
                if named and any(not _ISOTOPE.fullmatch(c) for c in channels):
                    raise InputError(
                        "Named MassHunter signal_columns must be isotope labels such as Fe56; "
                        "use an explicit generic CSV mapping for other named data columns."
                    )
                if qtegra and any(not _QTEGRA_ISOTOPE.fullmatch(c) for c in header[1:]):
                    raise InputError("Qtegra LDR requires mass-first isotope columns, such as 31P.")
                time_index = header.index(time_column)
                signal_indices = [header.index(c) for c in channels]
                awaiting_qtegra_metadata = qtegra
                continue

            if awaiting_qtegra_metadata:
                qtegra_metadata["channel_metadata"] = _qtegra_channel_metadata(row, header, row_number)
                awaiting_qtegra_metadata = False
                continue
            # Delimiter-only padding is accepted only immediately before a
            # recognized Printed footer. Interior/unterminated missing rows fail.
            if named and len(row) == len(header) and not any(cell.strip() for cell in row):
                if in_footer:
                    raise InputError(f"Row {row_number}: unexpected delimiter padding after MassHunter footer.")
                pending_delimited.append(row_number)
                pending_blank += 1
                continue
            blank = not row or (len(row) == 1 and not row[0].strip())
            if (named or qtegra) and blank:
                pending_blank += 1
                continue
            if named and row and row[0].strip().startswith("Printed:") and not any(cell.strip() for cell in row[1:]):
                if in_footer or not time_values:
                    raise InputError(f"Row {row_number}: unexpected MassHunter footer.")
                in_footer = True
                footer_padding_records.extend(pending_delimited)
                pending_delimited = []
                footer.append(row[0].strip())
                trailer_blanks += pending_blank
                pending_blank = 0
                continue
            if in_footer:
                raise InputError(f"Row {row_number}: nonblank content after MassHunter footer.")
            if pending_blank:
                raise InputError(f"Row {row_number}: blank rows interrupt the signal data.")
            if qtegra:
                if len(row) != len(header) + 1 or row[-1].strip():
                    raise InputError(f"Row {row_number}: Qtegra LDR requires a trailing empty field.")
                row = row[:-1]
            if len(row) != len(header):
                raise InputError(
                    f"Row {row_number}: expected {len(header)} fields, found {len(row)}."
                )
            value_raw = row[time_index]
            if time_kind == "numeric":
                value = _number(value_raw, time_column, row_number)
                scaled = value * _UNIT_SECONDS[time_unit]
                if not np.isfinite(scaled):
                    raise InputError(f"Row {row_number}: time overflows after conversion to seconds.")
            else:
                value = _datetime(value_raw, date_format, row_number)
                if first_value is not None and (
                        (value.utcoffset() is None) != (first_value.utcoffset() is None)):
                    raise InputError("Datetime values mix timezone-aware and naive timestamps.")
            if first_value is None:
                first_value, first_raw = value, value_raw
            if time_kind == "numeric":
                elapsed = (value - first_value) * _UNIT_SECONDS[time_unit]
            else:
                elapsed = (value - first_value).total_seconds()
            if not np.isfinite(elapsed):
                raise InputError(f"Row {row_number}: elapsed time is not finite.")
            if time_values and elapsed <= time_values[-1]:
                raise InputError(
                    f"Row {row_number}: time must be strictly increasing; "
                    "duplicate or out-of-order timestamps are not accepted."
                )
            values = [_number(row[index], name, row_number)
                      for index, name in zip(signal_indices, channels)]
            total = sum(values)
            if not np.isfinite(total):
                raise InputError(f"Row {row_number}: selected signal sum is not finite.")
            time_values.append(elapsed)
            signal_values.append(total)
    except csv.Error as exc:
        raise InputError(f"Malformed CSV near row {physical_offset + reader.line_num}: {exc}") from exc
    if header is None:
        raise InputError("Signal table header was not found.")
    if awaiting_qtegra_metadata:
        raise InputError("Qtegra LDR channel metadata row was not found.")
    if len(time_values) < 3:
        raise InputError("Signal data require at least three observations.")
    trailer_blanks += pending_blank
    if pending_delimited:
        raise InputError("Delimiter-only rows require a recognized trailing MassHunter Printed footer")
    metadata = {
        "recognized_footer_padding_physical_rows": footer_padding_records,
        "source": str(source),
        "sha256": digest,
        "format": source_format,
        "encoding": encoding,
        "delimiter": delimiter,
        "header_row": header_row,
        "time_column": time_column,
        "time_kind": time_kind,
        "time_unit": time_unit,
        "time_first": first_raw,
        "time_first_seconds": (first_value * _UNIT_SECONDS[time_unit]
                               if time_kind == "numeric" else None),
        "datetime_format": date_format,
        "output_time_unit": "s",
        "time_origin": "first_observation",
        "signal_columns": list(channels),
        "row_count": len(time_values),
        "preamble": preamble,
        "recognized_footer": footer,
        "trailing_blank_records": trailer_blanks,
    }
    if source_format == "csv":
        metadata["header_row_basis"] = "zero_based_physical_line"
        metadata["preamble_physical_text"] = preamble_physical_text
    if qtegra:
        metadata["qtegra"] = qtegra_metadata
        metadata["time_unit_source"] = "explicit_config"
        metadata["adapter_validation_scope"] = (
            "Five-point public developer fixture, Qtegra 2.7.2425.65, iCAP RQ, "
            "LDR File Version 1; not long-run or other-model validation.")
    return SignalData(np.asarray(time_values, dtype=float),
                      np.asarray(signal_values, dtype=float), metadata)
