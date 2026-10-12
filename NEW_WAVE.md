# New Wave event logs and short acquisitions

Use `laser.format: "new_wave_laser_log"` only for a New Wave event log with `Timestamp`, `Sequence Number`, `Comment`, and `Laser State` columns. In this format, a sequence number identifies a scan line or spot. Repeated `On` rows within that same sequence and sample update the state. One firing interval extends from the first `On` to the next `Off`. Crossing a sequence/sample boundary while firing, a missing `Off`, negative durations, or backwards firing-update timestamps remain input errors.

Choose an exact sample with `sample_selector`, or explicitly include the complete paired file with `sample_scope: "all"`. Do not set both. Sample selection can exclude intervening standards while retaining each selected interval's original times and source record numbers. Excluded records are listed as gaps; intervals are never joined across them. The source files are read-only.

An otolith configuration, saved next to the two original files, is:

```json
{
  "ms": {
    "path": "Otolith.csv",
    "format": "agilent_masshunter_csv",
    "signal_columns": ["Mg25", "P31", "Ca43", "Cu63", "Sr88", "Ba138", "Pb208", "Th232", "U238"],
    "signal_selection": "sum"
  },
  "laser": {
    "path": "Otolith_log.csv",
    "format": "new_wave_laser_log",
    "sample_selector": "Otolith"
  },
  "analysis": {"tolerance_s": 0.2, "coarse_window_s": 900},
  "export": {"mode": "rate_only"}
}
```

```text
python -B -X utf8 -m lagcheck run otolith_69.json --output otolith_69_results --export --preview Ca43
```

`"signal_selection": "sum"` fixes the matching signal to the sum of the listed columns. Without it, the default compares that sum with each listed column and uses the candidate with the highest starting-alignment peak.

For all 72 intervals, including three NIST612 standards, replace `sample_selector` with `"sample_scope": "all"` and use a new output directory. For the 96-spot example, select `Gabbro_data.csv`, `Gabbro_log.csv`, and `sample_scope: "all"`; explicitly list the desired matching channels from that CSV header and use `--preview Si29`. The recorded validation used the sum of all 9 otolith channels and all 24 gabbro channels for diagnosis. Choosing Ca43 or Si29 for display does not change the matching signal.

MassHunter comma-filled trailing rows are recognized as export padding only when immediately followed by a `Printed:` footer with empty remaining cells. Their physical row numbers are recorded. Interior missing rows, nonnumeric intensities, and unrecognized trailing content remain errors; no signal observations are silently deleted or replaced with zero.

## Conditional short-record analysis

New Wave inputs enable `analysis.short_acquisition` by default. Other input formats keep the previous standard behavior unless explicitly given `"short_acquisition": true`. Set this option to `false` to request the standard diagnostic alone.

The standard baseline has 10 blocks and requires at least 12 firing intervals per block, so it needs 120 intervals. The 15-block sensitivity setting needs 180 intervals. The conditional path uses at most the requested block count and at most `floor(interval_count / 12)` blocks. A rate still requires at least four usable blocks, correlation above 0.25, and peak margin above 0.01. Existing coverage, disagreement, held-out prediction and sensitivity thresholds are unchanged. A user-specified block count and changed gates are retained in the report, rather than silently overwritten.

If a short record has a nonboundary automatic correlation maximum above 0.25 but its uniqueness margin fails, that peak can be used for explicitly unverified conditional inspection. Its failed alignment check remains false. This does not establish the correct timing origin. A conditional rate, if available, remains provisional. This peak candidate is distinct from the existing retry at the input placement; no channel or signal is automatically changed.

`report.json` retains `standard_diagnostic`, the conditional checks and their results, and `short_acquisition.check_execution`. Each check is `pass`, `fail`, `incomplete`, or `not_executed`. Each sensitivity variant reports whether it was executed. Missing configurations are never treated as passes. A finite conditional rate has status `short_acquisition_conditional`; insufficient usable blocks retain a no-rate result. Both return exit code 2. Configuration, parsing and byte-integrity errors still return 1.

With `--export`, a finite conditional estimate can produce `provisional_laser_log.csv` and `provisional_manifest.json`. The timestamps are read back with the mapping held fixed, and the residual checks are saved. No short-path result is published as a verified corrected file while the full required sensitivity evidence is unavailable. Without a rate or explicit reference rate, no corrected log is manufactured. `--no-provisional` suppresses provisional log export and the existing input-placement retry; `--no-preview` suppresses plotting. `short_acquisition: false` disables the conditional short analysis, including its unverified peak inspection.

The laser-log screen also runs on short records and on spot or bracketed sequences. It can flag the intervals at which standards interrupt a sample sequence and irregular moves between spots. Such intervals belong to the acquisition sequence, and the screen cannot tell them from timestamp steps. Its state in `check_execution` is taken from `laser_log_screen` when no fit exists.

## Views and limits

In addition to `preview_pair.png`, the default New Wave preview writes:

- `new_wave_spatial.png`: actual first-On/Off coordinate endpoints. Moving lines use a stated straight-path, constant-progress display model between those endpoints. Spots use the interval mean at the recorded point, without spatial interpolation between spots. Separate sample labels remain separate panels for moving scans.
- `new_wave_signal_time.png`: the selected display channel over the complete MS time range and the mapped firing intervals.
- `new_wave_block_lags.png`: measured block lags and the baseline fit, when a fit exists.
- `new_wave_interval_signals.csv`: one row per interval and displayed mapping, including source record numbers, sample/sequence labels, sample counts and the chosen channel's interval mean.
- `new_wave_views.json`: source hashes, exact mapping source, coordinate assumptions and any geometry gaps. The right panel uses actual exported-log bytes when a log was written.

Repeated `On` records may contain stale coordinates, including a previous sample's last position at the same timestamp. They are retained in the audit but are not interpreted as measured trajectory vertices. A curved or otherwise unknown trajectory is not reconstructed from those updates. The coordinate preview is a transparent display model, not independent spatial accuracy or a tested iolite/HDIP import. The chosen plot channel is explicitly recorded; these views are not a multi-channel quantitative image export.

The 69-line sample-only example can still produce the original-coordinate image and diagnostics with no rate estimate. Such output means the input was read and inspected, not that a correction was obtained. The 72-interval and 96-spot examples can produce conditional estimates, images and provisional logs even though the full standard verification is unavailable. An external instrument's physical rate is not independently established by any of these outputs.
