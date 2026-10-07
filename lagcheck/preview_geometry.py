"""Shared channel loading and geometry helpers for visual inspection."""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from .models import InputError, IntegrityError
X_COLUMNS = ("X(um)", "Intended X(um)")
Y_COLUMNS = ("Y(um)", "Intended Y(um)")
SPEED_COLUMN = "Scan Velocity (um/s)"
DEFAULT_PITCH_UM = 15.0


def rasterize_true_y(t, chan, L, lags, x0, ncol, pix, y0, nrow, ypix):
    """Rows by stage y (not by line index): segments that share one y share one row.

    Each line is binned and gap-filled within its own sampled extent before it is written
    into the image, so the empty stretch between two segments of one row is never bridged.
    """
    img = np.full((nrow, ncol), np.nan)
    weights = np.zeros((nrow, ncol), dtype=np.int32)
    for d, lag in zip(L, lags):
        s, e = d["Start Time"] + lag, d["End Time"] + lag
        if e - s <= 0:
            continue
        i0, i1 = np.searchsorted(t, [s, e])
        if i1 - i0 < 1:
            continue
        frac = (t[i0:i1] - s) / (e - s)
        xpos = d["StartX"] + frac * (d["EndX"] - d["StartX"])
        col = np.clip(np.floor((xpos - x0) / pix).astype(int), 0, ncol - 1)
        acc = np.bincount(col, weights=chan[i0:i1], minlength=ncol)
        cnt = np.bincount(col, minlength=ncol)
        row = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
        v = np.flatnonzero(np.isfinite(row))
        if len(v) >= 2:
            inside = np.arange(v.min(), v.max() + 1)
            row[inside] = np.interp(inside, v, row[v])
        r = int(round((d["StartY"] - y0) / ypix))
        if 0 <= r < nrow:
            ok = np.isfinite(row)
            prior = weights[r, ok]
            img[r, ok] = (np.nan_to_num(img[r, ok])*prior+row[ok])/(prior+1)
            weights[r, ok] += 1
    return img


def _finite(value):
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def _sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8-sig"))
    except (OSError, ValueError) as exc:
        raise InputError(f"Cannot read {path}: {exc}") from exc


def _load_channel(config, report, channel):
    ms_config = {k: v for k, v in config["ms"].items() if k != "path"}
    path = Path(config["ms"]["path"])
    netcdf = ms_config.get("format") == "netcdf"
    key = "signal_variables" if netcdf else "signal_columns"
    stored = ms_config.get(key)
    if not isinstance(stored, list) or not stored:
        raise InputError(f"config.json ms.{key} must list the stored channels")
    if channel is None:
        channel = stored[0]
        note = ("same channel as the diagnosis" if len(stored) == 1 else
                f"the diagnosis matched the sum of {len(stored)} channels ({', '.join(stored)}); "
                "the preview shows the first stored channel only")
    else:
        if not isinstance(channel, str) or not channel.strip():
            raise InputError("Preview channel must be a column or variable name")
        note = ("same channel as the diagnosis" if stored == [channel] else
                "channel chosen for the preview; the diagnosis matched "
                + ("the sum of " if len(stored) > 1 else "") + ", ".join(stored))
    single = dict(ms_config, **{key: [channel]})
    if netcdf:
        from .netcdf import load_netcdf_signal
        ms = load_netcdf_signal(path, single)
    else:
        from .ms import load_signal
        ms = load_signal(path, single)
    if ms.metadata.get("sha256") != report["input_sha256"]["ms"]:
        raise IntegrityError("The MS input no longer matches the hash recorded in report.json; no preview")
    return channel, stored, note, ms


def _speeds_from_log(table, header):
    """Scan speed at each firing start, inheriting the last nonblank cell (ESL convention)."""
    if SPEED_COLUMN not in header:
        return None
    col = header.index(SPEED_COLUMN)
    starts = {on for on, _ in table["interval_rows"]}
    current, speeds = None, {}
    for i, row in enumerate(table["rows"]):
        cell = row[col].strip() if col < len(row) else ""
        if cell:
            try:
                current = float(cell)
            except ValueError:
                current = None
        if i in starts:
            speeds[i] = current
    return [speeds[on] for on, _ in table["interval_rows"]]


def _pitch(positions):
    unique = np.unique(np.round(np.asarray(positions, float), 6))
    gaps = np.diff(unique)
    gaps = gaps[gaps > 0]
    if gaps.size:
        return float(np.median(gaps)), "median spacing of the distinct line positions in the log"
    return DEFAULT_PITCH_UM, "default 15 um (only one line position in the log)"


def _geometry(laser, ms, settings):
    """Map this package's lines to StartX/EndX/StartY for the rasterizer."""
    table, lines = laser.table, laser.lines
    header = [h.strip() for h in table["header"]]
    rows, pairs = table["rows"], table["interval_rows"]
    n = len(lines)
    log_speeds = _speeds_from_log(table, header)
    declared = _finite((settings or {}).get("scan_speed_um_s"))
    speed, speed_source = None, "not available (no analysis.scan_speed_um_s and no usable log column)"
    if declared is not None and declared > 0:
        speed, speed_source = declared, "analysis.scan_speed_um_s"
    elif log_speeds and log_speeds[0] is not None and log_speeds[0] > 0:
        speed, speed_source = log_speeds[0], f"laser log column {SPEED_COLUMN!r}, first firing line"
    uniform = None if not log_speeds else len(set(log_speeds)) == 1
    if laser.metadata.get('format') == 'new_wave_laser_log':
        speed, uniform = None, None
        speed_source = ('New Wave metadata can inherit positioning speed at the first On; '
                        'no verified firing speed or spot-speed displacement conversion is available')
    xcol = next((c for c in X_COLUMNS if c in header), None)
    ycol = next((c for c in Y_COLUMNS if c in header), None)
    note = None
    if xcol and ycol and laser.metadata.get("format") in ("laser_state_csv", "new_wave_laser_log"):
        ix, iy = header.index(xcol), header.index(ycol)
        try:
            sx = np.array([float(rows[a][ix]) for a, _ in pairs])
            ex = np.array([float(rows[b][ix]) for _, b in pairs])
            sy = np.array([float(rows[a][iy]) for a, _ in pairs])
            ey = np.array([float(rows[b][iy]) for _, b in pairs])
        except (ValueError, IndexError):
            note = f"log columns {xcol!r}/{ycol!r} exist but are not numeric at every firing endpoint"
        else:
            if not all(np.isfinite(v).all() for v in (sx, ex, sy, ey)):
                note = f"log columns {xcol!r}/{ycol!r} contain non-finite positions"
            elif np.abs(ex - sx).max() == 0 and np.abs(ey - sy).max() == 0:
                note = f"log columns {xcol!r}/{ycol!r} show no movement within any firing interval"
            else:
                swapped = bool(np.median(np.abs(ey - sy)) > np.median(np.abs(ex - sx)))
                along_start, along_end, across = (sy, ey, sx) if swapped else (sx, ex, sy)
                along_change = np.abs(along_end-along_start)
                across_change = np.abs(ex-sx) if swapped else np.abs(ey-sy)
                # Endpoint coordinates rounded to micrometres may differ slightly.
                # A mixed/diagonal trajectory is not a horizontal raster line.
                tolerance = np.maximum(1e-6, along_change*1e-6)
                if np.any(across_change > tolerance):
                    note = 'diagonal or mixed-direction firing intervals: line-time view used; endpoints do not define a common raster axis'
                else:
                    pitch, pitch_source = _pitch(across)
                    return dict(source=f"laser log columns {xcol!r} and {ycol!r} at the firing On/Off records",
                                swapped=swapped, unit="um", along_start=along_start, along_end=along_end,
                                across=across, pitch=pitch, row_pitch=pitch, pitch_source=pitch_source,
                                note=None), speed, speed_source, uniform
    elif laser.metadata.get("format") in ("laser_state_csv", "new_wave_laser_log"):
        note = f"no {X_COLUMNS[0]}/{Y_COLUMNS[0]} columns in the log"
    else:
        note = "interval_csv input carries no stage positions"
    durations = np.array([float(l["End Time"] - l["Start Time"]) for l in lines])
    if speed is not None:
        geometry = dict(source="synthetic rows: x = scan speed * time along the line, one line per row",
                        swapped=False, unit="um", along_start=np.zeros(n), along_end=speed * durations,
                        across=np.arange(n, dtype=float) * DEFAULT_PITCH_UM, pitch=DEFAULT_PITCH_UM,
                        row_pitch=DEFAULT_PITCH_UM, pitch_source="default 15 um", note=note)
    else:
        t = np.asarray(ms.time, float)
        dt = float((t[-1] - t[0]) / (len(t) - 1))
        geometry = dict(source="synthetic rows: x = time along the line in seconds, one line per row",
                        swapped=False, unit="s", along_start=np.zeros(n), along_end=durations,
                        across=np.arange(n, dtype=float), pitch=dt, row_pitch=1.0,
                        pitch_source="mean MS sample spacing (one sample per pixel)", note=note)
    return geometry, speed, speed_source, uniform
