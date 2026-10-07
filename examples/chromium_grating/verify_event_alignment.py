"""Recalculate the documented first-Si29-event offset without tuning it.

Exit 0: matched; exit 2: mismatch; exit 1: input or tool error.
"""
from pathlib import Path
import json
import math
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))


def main():
    try:
        import numpy as np
        from lagcheck.ms import load_signal
        from lagcheck.laser import load_laser
        config = json.loads((HERE / "first_si_event.json").read_text(encoding="utf-8"))
        reference = json.loads((HERE / "coordinate_reference.json").read_text(encoding="utf-8"))
        ms = load_signal(HERE / config["ms"]["path"], {k: v for k, v in config["ms"].items() if k != "path"})
        laser = load_laser(HERE / config["laser"]["path"], {k: v for k, v in config["laser"].items() if k != "path"})
        raw_time = ms.time + ms.metadata["time_first_seconds"]
        floor = float(np.median(ms.signal[raw_time < 10.0]))
        plateau = float(np.median(ms.signal[(raw_time > 25.0) & (raw_time < 40.0)]))
        threshold = float(np.sqrt(floor * plateau))
        rising = np.flatnonzero((raw_time[:-1] < 45.0) & (raw_time[1:] < 45.0)
                               & (ms.signal[:-1] < threshold) & (ms.signal[1:] >= threshold))
        if not len(rising):
            print(json.dumps({"status": "mismatch", "reason": "No event in the declared search window"}))
            return 2
        j = int(rising[0])
        fraction = float((threshold - ms.signal[j]) / (ms.signal[j+1] - ms.signal[j]))
        crossing = float(ms.time[j] + fraction * (ms.time[j+1] - ms.time[j]))
        first_on = float(laser.lines[0]["Start Time"])
        offset = crossing - first_on
        values = {"ms_first_sample_seconds": ms.metadata["time_first_seconds"],
                  "laser_first_on_native_s": first_on, "n_intervals": len(laser.lines),
                  "n_signal_points": len(ms.time), "floor_median_cps": floor,
                  "plateau_median_cps": plateau, "geometric_mean_threshold_cps": threshold,
                  "first_si_crossing_native_s": crossing,
                  "first_si_crossing_raw_s": crossing + ms.metadata["time_first_seconds"],
                  "alignment_offset_s": offset}
        differences = [key for key, value in values.items()
                       if not math.isclose(value, reference[key], rel_tol=0.0, abs_tol=1e-9)]
        if not math.isclose(offset, config["analysis"]["alignment_offset_s"], rel_tol=0.0, abs_tol=1e-9):
            differences.append("configured alignment_offset_s")
        print(json.dumps({"status": "matched" if not differences else "mismatch",
                          "differences": differences, "values": values}, indent=2))
        return 0 if not differences else 2
    except Exception as exc:
        print(json.dumps({"status": "tool_or_input_error", "error": type(exc).__name__ + ": " + str(exc)}))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
