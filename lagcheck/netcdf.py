"""Explicit variable mapping for classic NetCDF signals, independent of vendor.

This reader has synthetic-file tests only. It does not establish compatibility
with any particular PerkinElmer, Syngistix or NexION export. Variable names and
time units must be supplied after inspecting the actual exported file.
"""
from __future__ import annotations

from collections.abc import Mapping
import hashlib
from pathlib import Path
import struct

import numpy as np

from .models import InputError, SignalData


_UNITS = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "min": 60.0}
_DEFAULT_FILL = {("i", 1): -127, ("i", 2): -32767, ("i", 4): -2147483647,
                 ("f", 4): 9.969209968386869e36, ("f", 8): 9.969209968386869e36}


def _open(path):
    try:
        from scipy.io import netcdf_file
    except ImportError as exc:
        raise InputError("Classic NetCDF input requires scipy") from exc
    try:
        source = Path(path).resolve()
        with source.open("rb") as stream:
            signature = stream.read(4)
        if signature not in (b"CDF\x01", b"CDF\x02"):
            raise InputError("Only classic NetCDF CDF1/CDF2 files are supported; NetCDF4/HDF5 is not supported")
        digest = _hash(source)
        return source, netcdf_file(source, mode="r", mmap=False, maskandscale=False), digest
    except InputError:
        raise
    except (OSError, TypeError, ValueError, OverflowError, KeyError, IndexError, EOFError, struct.error) as exc:
        raise InputError(f"Cannot open classic NetCDF file {path!s}: {exc}") from exc


def _hash(path):
    digest = hashlib.sha256()
    try:
        with Path(path).open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise InputError(f"Cannot hash NetCDF input: {exc}") from exc
    return digest.hexdigest()


def _summary(value):
    """JSON-safe bounded metadata; never return an unbounded attribute array."""
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        return value if len(value) <= 512 else {"preview": value[:512], "length": len(value), "truncated": True}
    if isinstance(value, np.ndarray):
        if value.size > 16:
            return {"shape": list(value.shape), "dtype": str(value.dtype),
                    "preview": [_summary(v) for v in value.flat[:16]], "truncated": True}
        return [_summary(v) for v in value.flat]
    if isinstance(value, np.generic):
        return _summary(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    if value is None or isinstance(value, (int, float, bool)):
        return value
    return str(value)[:512]


def inspect_netcdf(path) -> dict:
    """List shapes, dtypes, units and bounded global metadata without values."""
    source, dataset, digest = _open(path)
    with dataset:
        variables = {name: {"shape": list(variable.shape), "dimensions": list(variable.dimensions),
                            "dtype": str(variable.data.dtype),
                            "units": _summary(variable._attributes.get("units")),
                            "attribute_names": sorted(variable._attributes)}
                     for name, variable in dataset.variables.items()}
        result = {"source": str(source), "sha256": digest, "format": "classic_netcdf",
                  "dimensions": dict(dataset.dimensions), "variables": variables,
                  "global_attributes": {str(k): _summary(v) for k, v in dataset._attributes.items()},
                  "vendor_validation": "generic mapping only; no NexION export has been validated"}
    if _hash(source) != digest:
        raise InputError("NetCDF input changed while inspecting")
    result["input_unchanged"] = True
    return result


def _array(dataset, name):
    if name not in dataset.variables:
        raise InputError(f"NetCDF variable {name!r} does not exist")
    variable = dataset.variables[name]
    attributes = variable._attributes
    if "scale_factor" in attributes or "add_offset" in attributes:
        raise InputError(f"Variable {name!r} uses packed scale_factor/add_offset metadata; unpack explicitly before import")
    if "_Unsigned" in attributes:
        raise InputError(f"Variable {name!r} uses an unsupported _Unsigned convention")
    raw = variable.data
    if np.ma.isMaskedArray(raw) and np.any(np.ma.getmaskarray(raw)):
        raise InputError(f"Variable {name!r} contains masked values")
    raw = np.asarray(raw)
    if raw.ndim != 1:
        raise InputError(f"Variable {name!r} must be one-dimensional, found shape {raw.shape}")
    if raw.dtype.kind not in "iuf":
        raise InputError(f"Variable {name!r} must contain real numeric values")
    if not np.all(np.isfinite(raw)):
        raise InputError(f"Variable {name!r} contains NaN or infinite values")
    sentinels = []
    for key in ("_FillValue", "missing_value"):
        if key in attributes:
            supplied = np.asarray(attributes[key])
            if supplied.dtype.kind not in "iuf" or supplied.size == 0:
                raise InputError(f"Variable {name!r} has unsupported {key} metadata")
            sentinels.extend(supplied.ravel())
    default = _DEFAULT_FILL.get((raw.dtype.kind, raw.dtype.itemsize))
    if default is not None:
        sentinels.append(np.asarray(default, dtype=raw.dtype).item())
    if any(np.any(raw == value) for value in sentinels):
        raise InputError(f"Variable {name!r} contains a fill or missing value")
    result = raw.astype(np.float64, copy=True)
    if not np.all(np.isfinite(result)):
        raise InputError(f"Variable {name!r} overflows floating-point conversion")
    return result


def load_netcdf_signal(path, config) -> SignalData:
    """Load explicitly mapped 1-D channels, sum them and return relative seconds.

    Required configuration keys are ``time_variable``, ``time_unit`` (s, ms, us,
    min), and a nonempty list ``signal_variables``. No variables or units are
    inferred from vendor names or metadata. Arrays are never sorted, flattened,
    interpolated, stripped of missing values or silently unpacked.
    """
    if not isinstance(config, Mapping):
        raise InputError("NetCDF configuration must be a mapping")
    for key in ("time_variable", "time_unit", "signal_variables"):
        if key not in config:
            raise InputError(f"NetCDF configuration requires {key!r}")
    time_name, unit, names = config["time_variable"], config["time_unit"], config["signal_variables"]
    if not isinstance(time_name, str) or not time_name:
        raise InputError("time_variable must be a nonempty variable name")
    if not isinstance(unit, str) or unit not in _UNITS:
        raise InputError("time_unit must be explicitly set to s, ms, us or min")
    if (not isinstance(names, list) or not names
            or any(not isinstance(name, str) or not name for name in names)
            or len(names) != len(set(names))):
        raise InputError("signal_variables must be a nonempty list of distinct variable names")
    if time_name in names:
        raise InputError("time_variable cannot also be a signal variable")
    source, dataset, digest = _open(path)
    with dataset:
        if time_name not in dataset.variables:
            raise InputError(f"NetCDF variable {time_name!r} does not exist")
        time_variable=dataset.variables[time_name]
        declared=time_variable._attributes.get('units')
        if isinstance(declared,bytes):
            declared=declared.decode('ascii',errors='replace')
        aliases={'s':'s','sec':'s','second':'s','seconds':'s','ms':'ms','millisecond':'ms',
                 'milliseconds':'ms','us':'us','microsecond':'us','microseconds':'us',
                 'min':'min','minute':'min','minutes':'min'}
        known_unit=aliases.get(declared.lower().strip().split(' since ')[0]) if isinstance(declared,str) else None
        if known_unit is not None and known_unit!=unit:
            raise InputError(f"Configured time_unit {unit!r} conflicts with NetCDF time units {declared!r}")
        time = _array(dataset, time_name)
        if len(time) < 3:
            raise InputError("NetCDF signal requires at least three observations")
        if not np.all(time[1:] > time[:-1]):
            raise InputError("NetCDF time must be strictly increasing")
        arrays = [_array(dataset, name) for name in names]
        dimensions={name:list(dataset.variables[name].dimensions) for name in [time_name]+names}
        if any(dataset.variables[name].dimensions!=time_variable.dimensions for name in names):
            raise InputError("NetCDF time and signal variables must share the same named dimension")
        if any(len(array) != len(time) for array in arrays):
            raise InputError("NetCDF time and selected signal variables must have equal lengths")
        with np.errstate(over="ignore", invalid="ignore"):
            elapsed = (time - time[0]) * _UNITS[unit]
            first_seconds = time[0] * _UNITS[unit]
            signal = np.sum(np.stack(arrays), axis=0)
        if not np.all(np.isfinite(elapsed)) or not np.isfinite(first_seconds):
            raise InputError("NetCDF time overflows conversion to seconds")
        if not np.all(elapsed[1:] > elapsed[:-1]):
            raise InputError("NetCDF time precision is insufficient after conversion to seconds")
        if not np.all(np.isfinite(signal)):
            raise InputError("Selected NetCDF signal sum is not finite")
        units = {name: _summary(dataset.variables[name]._attributes.get("units"))
                 for name in [time_name] + names}
    if _hash(source) != digest:
        raise InputError("NetCDF input changed while reading")
    metadata = {"source": str(source), "sha256": digest, "format": "classic_netcdf",
                "time_variable": time_name, "time_kind": "numeric", "time_unit": unit,
                "time_first": float(time[0]), "time_first_seconds": float(first_seconds),
                "output_time_unit": "s", "time_origin": "first_observation",
                "signal_variables": list(names), "signal_columns": list(names),
                "signal_operation": "sum of explicitly selected variables",
                "row_count": len(time), "variable_units_metadata": units,
                "variable_dimensions": dimensions,
                "time_unit_metadata_checked": known_unit is not None,
                "input_unchanged": True,
                "vendor_validation": "generic mapping only; no NexION export has been validated"}
    return SignalData(elapsed, signal, metadata)
