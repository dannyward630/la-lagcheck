# Troubleshooting when no rate is available or a check is not met

Companion guide for LA-LagCheck 0.8.0. Run the commands from the package directory with its requirements installed. Keep configuration copies beside the original configuration so relative input paths retain their meaning. The guide covers four situations: no rate is estimated, a timing check is not met, no provisional log is written from the estimate, and a tolerance in micrometres is refused.

Open `report.json` and `preview_pair.json` first. Read `reason`, `alignment`, `rate_ppm`, `warnings`, and the mapping source. A raw or unregistered view gives you something to inspect; it does not mean a correction was estimated. If `INCOMPLETE.json` exists, the run has an input or execution error. Its raw inspection plot does not make a corrected log usable. For each check that is not met, `report.md` gives the measured value and its limit where the report holds them. `notes` in `report.json`, repeated in `report.md` and in the printed line, states when the record has too few firing intervals, when the iteration did not converge and whether a supplied reference rate was used.

Use a new configuration copy and a new output directory for each trial. Keep the original input files and previous results. The program first retains its strict diagnostic. With provisional output enabled, an automatic input-placement retry may supply a provisional rate under the conditions below. The later configuration edits are explicit user choices.

## Understand the automatic input-placement retry

If you have not declared a starting alignment, the program can repeat diagnosis when automatic alignment is unidentifiable or its rate fails the agreement or held-out prediction checks. It uses the parent input placement (`alignment_offset_s = 0` for an unsegmented input) with the same matching signal and other settings. Segment offsets preserve that placement by accounting for the shift from the parent origin; each segment is not reset to the start of the mass-spectrometer record. It does not choose a different channel. The retry is used only if the candidate rate is finite, the time-scale factor is positive, and agreement and prediction explicitly pass at both the initial-fit and iteration stages. Missing checks do not pass. Iteration convergence is reported but is not an additional requirement for provisional adoption. Its correction is provisional, and the strict result remains in the report.

Read the top-level `status`, `reason`, `rate_ppm`, `alignment`, `checks` and `input_diagnostic_status` as the unchanged strict result. Read `alignment_fallback` separately: it records `attempted`, `adopted`, `status`, `triggered_reason`, `alignment_offset_s`, `required_checks`, `failed_adoption_checks`, `diagnostic` and `note`. The candidate diagnostic is also saved in `alignment_fallback_report.json` and its block files. Adoption is marked `provisional_adopted` with rate source `input_position_fallback`; the run returns exit code 2. Failed agreement or prediction checks remain visible.

To disable this automatic retry while retaining other provisional outputs, run:

```text
python -B -X utf8 -m lagcheck run config.json --output ../trial_without_fallback --export --no-alignment-fallback
```

The same switch is available in `wizard`. `--no-provisional` also disables the retry; `--no-preview` only disables plots. A command without `--export` may still run the retry and show its provisional mapping, but writes no corrected log. An adopted retry takes precedence over the original estimate and an explicitly supplied reference. The rate of an adopted retry is written only when it is eligible for provisional export (see "When no provisional log is written from the estimate" below); when it is not, a supplied reference rate is written instead. If the retry is rejected, an original rate can still be shown or exported under the existing provisional rules, with both the original failed checks and the rejected retry checks recorded.

Zero refers to the placement supplied by the readers or saved workspace, not a measurement that the instruments started together. With the default ±25 s block-search window, the retry assumes corresponding events fall within that window. Neither simultaneous-looking record starts nor passing agreement and prediction checks proves the physical starting alignment. Use an identifiable starting event or independent marker to assess absolute alignment.

If the retry is rejected and no other eligible rate is available, the program retains the original/declared-alignment view and records the failed retry checks. Input and time-grid errors are not bypassed. The manual options below remain available when justified by the data.

## Declare a starting alignment when you have a basis for it

If the automatic starting alignment is not identifiable but you have established the matching start independently, copy the existing configuration to `config_alignment.json`. In its existing `analysis` object, set `alignment_offset_s` to that value in seconds. Keep all other required configuration entries.

```json
"analysis": {
  "tolerance_s": 0.2,
  "alignment_offset_s": 3.0
}
```

The value `3.0` is a syntax example, not a suggested offset. The program places a laser time `t` at MS time `t + alignment_offset_s`, before estimating a rate. Both readers rebase their coordinates to their own first retained timestamps; the offset must refer to those reader coordinates. It is not automatically the difference between acquisition headers. Record where the value came from with the trial. Do not declare zero merely to make the analysis continue.

```text
python -B -X utf8 -m lagcheck run config_alignment.json --output ../trial_alignment --export
```

This declaration sets the initial diagnostic alignment. It does not by itself establish the first MS sample's absolute time for an affine export or verify the final reconstruction.

### What counts as an independent basis

A basis is independent when it is obtained without the alignment choice that it is meant to justify. It does not have to be a reference material or a second instrument. These can serve:

- A first firing that can be identified without doubt after recorded background, paired with the one signal response that belongs to it. Allow for the transport and response delay between the firing and the signal.
- A marker that is recorded on both time axes.
- A feature of a separately registered image or of the stage coordinates, when the feature does not repeat and its correspondence to the measured signal is documented. Establish the correspondence apart from the reconstructed image that is being assessed. A boundary that only looks the same is not enough when its signal response cannot be identified uniquely.

Record the original event times, the conversion from the time origin of each reader, the part of the signal response that was used (for example its onset) and the uncertainty in locating it. Explain how the pair of events was identified. Do not select the offset that best fits the repeated scan pattern. A hardware marker relates the two clocks. It does not say at which firing position a measured signal belongs.

These do not count as an independent basis:

- an offset saved in a workspace;
- two records that start at the same time, or an offset of zero chosen for that reason;
- the highest correlation peak of a repeated scan pattern;
- the first point of the same fitted lag line, extrapolated back to the first line;
- a reconstruction that looks better;
- a similar average rate obtained from the same raw data with another algorithm. It does not establish the starting position or the geometry of the whole scan.

## Read the laser-log screen

`laser_log_screen.flagged` lists line-change intervals that differ from the median of their neighbours by more than the limit. The check `no_flagged_line_change_interval` is then not met and the result is not supported, whatever the other checks show. A flag marks a line-change interval that is unusually long or short relative to its neighbours. It may reflect an acquisition pause, a timestamp step or variation in the acquisition sequence, and the log alone does not establish the cause. Compare the signal on either side. A change in lag close to the interval excess supports a timestamp step; similar lags are consistent with a real change in the line-change interval. Interpret either result together with the acquisition sequence and the signal-matching uncertainty. The program does not divide a record at a flagged interval or remove the excess. `python -B -X utf8 -m lagcheck screen-log config.json` screens a log without a signal file. A known pause can be accepted by raising `analysis.line_change_limit_s` for that run and stating the reason with the result. The raised limit and the number of intervals that the default limit would flag are written to the report and to the manifest of any exported log. Intervals that recur by design, such as bracketing standards, are flagged as well. Groups of similarly displaced adjacent intervals can escape detection because they shift the local reference. Random scatter of the line-change intervals can also produce a flag once its standard deviation reaches a few tenths of a second: in 200 simulated logs of 600 independent intervals for each value, with normally distributed scatter about a 3 s mean, no interval was flagged at standard deviations of 0.05, 0.1 or 0.2 s, and 6 and 44 logs had a flagged interval at 0.3 s and 0.5 s.

## Blocks at the start or end of the signal record

When the signal record begins as the signal of the first firing arrives, the lag search of the first block cannot extend before the record. Once the rate has been removed in the iteration, the first-block peak can lie at that boundary, which fails `peak_not_at_search_boundary` or `closure_peak_not_at_boundary`, or, when that block is already clipped by the record, `all_sensitivity_evidence_supported` in the finer block divisions. The program now handles this by default (`analysis.skip_record_edge_lines`, default `true`). If a peak lies at a boundary imposed by the start or end of the record, it repeats the diagnosis without the scan lines that start within 30 s of the first recorded sample, or end within 30 s of the last, and leaves the same lines out of every fit. Open `record_edge_lines` in the report: `applied` says whether the repeat was made, `sides`, `skipped_at_start` and `skipped_at_end` say which lines were left out, and `guard_s` gives the 30 s. Coverage is still measured against all lines. The rule uses the neighbouring lines and the laser-off intervals between lines as the missing signal, so it assumes that the record contains them.

If a peak is at a record boundary but `record_edge_lines.applied` is `false`, the rule could not be applied: no line lies within the guard of that end, or fewer than twelve lines per block would remain (`not_applied_reason`, `available_at_start` and `available_at_end` in the report). The check then stays unmet. `record_boundary_blocks` in the report names the blocks concerned. Recording background before the first firing, and after the last, avoids the condition. For records already acquired, `"exclude_record_boundary_blocks": true` under `analysis` is a declared alternative: such a block leaves the fit but must still agree with the fitted line, the status becomes `supported_excluding_record_boundary_blocks`, and an exported log stays provisional. The agreement test of an excluded block cannot see a lag beyond the record boundary, so use the option only when the record is known to contain the signal of the first firing. Setting `"skip_record_edge_lines": false` restores the 0.6.0 behaviour, for example to reproduce an earlier result.

## Explicitly choose a different matching signal

For the CSV signal formats the program chooses the matching signal itself (`ms.signal_selection`, default `"highest_alignment_peak"`): it compares the sum of the configured columns and each configured column by the peak of the automatic starting-alignment search, keeps the sum unless a single column is higher by more than 1e-9, and records every candidate under `signal_selection` in the report. A candidate can win on peak height without meeting the starting-alignment criteria, so read that record. NetCDF input keeps `"sum"`. To use the sum without comparison, set `"signal_selection": "sum"` under `ms`. To force one channel, copy the configuration to `config_single_channel.json` and change the existing CSV `ms.signal_columns` list. For example, if the input contains an exact `Fe56` column:

```json
"signal_columns": ["Fe56"]
```

For a NetCDF mapping, change its existing `signal_variables` list instead, using an actual variable name from that file. Keep its time mapping and acquisition selection unchanged. With `"signal_selection": "sum"` a list of multiple channels is summed for diagnosis; a list of one channel makes that channel the matching signal.

```text
python -B -X utf8 -m lagcheck run config_single_channel.json --output ../trial_single_channel --export
```

`--preview Fe56` changes only the displayed channel. It does not change the signal used to estimate the rate. Retain the channel choices and all trial reports; inspect why the diagnostic changed rather than selecting only a passing result.

## Supply a measured reference rate when the acquisition has no estimate or the estimate is withheld

If a rate measured for the same instrument pair is appropriate and you elect to use it, supply it explicitly:

```text
python -B -X utf8 -m lagcheck run config.json --output ../trial_reference --export --reference-rate-ppm 170.9
```

`170.9` is the study's instrument-pair example, not a universal default. Replace it with your justified reference value and retain its source. A reference rate must come from a documented measurement on the same instrument pair, under conditions that justify using it for this acquisition. An adopted input-placement retry takes precedence over both the original strict estimate and a supplied reference. The reference is used when no acquisition estimate is available after this selection, or when that estimate is not eligible for provisional export (see "When no provisional log is written from the estimate" below). It does not replace an eligible estimate, supply a missing time origin, establish the starting alignment, repair invalid input, or turn an unsupported diagnostic into a supported one. It is written as a rate-only mapping. It takes no fitted constant from a withheld estimate, also when an affine export was requested. When a reference rate is supplied, `reference_rate` in `report.json` records whether it was used and the reason (`no_acquisition_estimate`, `acquisition_estimate_not_eligible`, `acquisition_estimate_available`, `export_not_requested` or `provisional_output_disabled`), and `report.md` and the printed line repeat that statement.

When input and export requirements are satisfied, this route can produce `provisional_laser_log.csv` and `provisional_manifest.json`, marked not verified with the reference rate source. Without a valid starting alignment, fixed-mapping read-back can remain `not_run`; the display says that its clock origins are unregistered. Do not treat the plot as confirmation of absolute alignment.

## When a check is not met: what a provisional log is for

A timing check that is not met remains a limitation of the result when a provisional log is available. A provisional log can be used to inspect a candidate mapping and to compare it with the original reconstruction. Keep its source and its failed checks with every use. It does not support a statement about positions by itself. For such a statement, compare the result with independent spatial evidence that meets the accuracy you need, over the part of the scan that you use. Agreement at the start of the scan does not establish agreement at the end. Without such evidence, keep the result as a candidate that is not verified and make no quantitative claim about its spatial accuracy. A preview that was drawn is not a check of spatial accuracy.

Before you use a provisional log, read the entries below in `report.json`. `report.md` gives a summary of them, and the JSON files hold the complete records.

Start with the rate that the file applies and its source: `provisional_rate_ppm` and `provisional_rate_source`. The source says where items 2 to 7 are read:

- `acquisition_estimate`: read them at the top level of `report.json`.
- `input_position_fallback`: the log applies the rate of the adopted input-placement retry. Read the fields of the same names in `alignment_fallback.diagnostic`, for example `alignment_fallback.diagnostic.baseline.rate_ppm`. `alignment_fallback_report.json` holds the same diagnostic. The top level of `report.json` is then the strict diagnostic. It describes another estimate and is not evidence for this log.
- `user_reference_rate`: the log applies the reference rate that you supplied, not an estimate of this acquisition. Items 2 to 7 then describe the estimate that was not written. They say why it was not used, and they are not evidence for the reference rate. `provisional_export_eligibility.candidate.source` names that estimate: read it at the top level for `acquisition_estimate` and in `alignment_fallback.diagnostic` for `input_position_fallback`. When `provisional_export_eligibility.status` is `no_acquisition_rate`, there is no estimate to read. Keep the source of the reference rate and its uncertainty with the result.

1. The rate that the file applies and its source, as above. `provisional_export_eligibility.applied_mapping` holds the same record, and the preview reads it.
2. The single-pass rate and the standard error of its slope: `baseline.rate_ppm` and `baseline.slope_se_ppm`. The standard error belongs to the single-pass fit. It is not a standard error of the rate after iteration.
3. The number of retained blocks and the coverage: `baseline.n_retained` and `coverage_fraction`.
4. The largest held-out prediction error: `baseline.heldout_max_s`.
5. The blocks whose lag differs from the fitted line by more than the tolerance: `contradictory_high_quality_blocks`.
6. The rates under the alternative settings, their range and the settings that gave no rate: `sensitivity`, `sensitivity_envelope_s` and `provisional_export_eligibility.missing_sensitivity_configs`. The sign of a missing estimate is unknown.
7. Whether the iteration converged: `closure.converged`. When it did not and the log applies an estimate of this acquisition, that estimate is the single-pass rate. The checks named `closure_*` were then evaluated at the last trial rate of the iteration and are not evidence for the single-pass rate. A log that applies a reference rate applies that rate, whether the iteration converged or not.
8. The laser-log screen: `laser_log_screen.flagged`, at the top level of `report.json`.
9. The read-back of the written file with its mapping held fixed. Where it is recorded depends on how the provisional log was written:
   - `provisional_failure_stage` is `readback` and `provisional_rate_source` is `acquisition_estimate`: the provisional log is the checked candidate, byte for byte. Its read-back is `post_export_verification` in `report.json`, and the diagnostic of that read-back is `post_export_report.json`. `provisional_readback` repeats the status and the failed checks of that read-back. No `provisional_readback_report.json` is written.
   - Otherwise the provisional log was written as its own file. Its read-back is `provisional_readback` in `report.json`. When that read-back was run, its diagnostic is `provisional_readback_report.json`. This includes a reference-rate log that was written after the read-back of the candidate failed. `post_export_verification` is then the read-back of that candidate, not of the provisional log.
   - When `provisional_readback.status` is `not_run`, no starting alignment was available and no read-back was made. `readback_status` in `provisional_manifest.json` is then `not_run` as well, and the run directory holds no `provisional_readback_report.json`.

   Check that the record refers to the file and the rate that you use. `provisional_manifest.json` gives the hash of the file (`output_sha256`), its rate (`rate_ppm`) and the source of the rate (`rate_source`). `readback_record` in `preview_pair.json` names the record that the preview read.
10. File integrity and format: the run directory holds no `INCOMPLETE.json`, and the hash in the manifest is that of the file.

Compare each time difference with the tolerance that you chose for the intended use. Do not rely on the fitted rate alone. A residual rate becomes a time difference when its magnitude in µs s⁻¹ is multiplied by the firing span in seconds and by 10⁻⁶. The result is the lag change that belongs to that residual slope. It is not a bound on the total alignment error. The firing span is the time from the first firing start to the last firing end (`duration_s`).

The internal read-back does not test the receiving software. Inspect the imported timing and coordinates there as well. Keep the original data, the settings and the diagnostic report with any reconstruction made from a provisional file.

## When no provisional log is written from the estimate

The rate estimated from the acquisition is written to a provisional log only when it is eligible for provisional export. `provisional_export_eligibility` in `report.json` records this for the diagnostic that an export uses: the adopted input-placement retry, otherwise the main diagnostic. The eligibility restricts provisional output only. The rate, the checks, the status and the exit code are unchanged, and a supported result is still exported through the verified path.

The estimate is withheld when one of three conditions holds. T is the firing span in seconds, and the tolerance is that of the run in seconds.

1. `rate_direction`: the direction of the single-pass rate is uncertain. The single-pass rate lies within its screening half-width of zero, and the half-width times T times 10⁻⁶ exceeds the tolerance. The half-width is the 0.975 quantile of Student's t with n - 2 degrees of freedom times the standard error of the single-pass slope, where n is the number of retained blocks. It is a scale of the conditional fit, not a confidence interval of the reported rate.
2. `sign_across_settings`: the single-pass rate, the candidate rate and the rates of the alternative settings differ in sign, and their range times T times 10⁻⁶ exceeds the tolerance.
3. `unconverged_iteration`: the iteration did not converge, and the largest difference between the rate to be written and a trial rate of the iteration, or the residual rate at the last trial, times T times 10⁻⁶ exceeds the tolerance.

A condition is triggered only above the tolerance. An effect equal to the tolerance does not trigger it. A negative rate is not a reason to withhold an estimate, and a rate near zero is not withheld when its uncertainty has no effect above the tolerance.

`status` has four values:

- `not_blocked`: none of the three conditions was triggered. This says nothing more about the estimate. The failed checks and the missing estimates in the report still apply.
- `withheld`: at least one condition was triggered. `triggered_conditions` and `reasons` name it with its numbers.
- `not_assessable`: no condition was triggered, but a quantity that a condition needs is missing or not finite. The estimate is not written. A missing quantity is never read as zero.
- `no_acquisition_rate`: no rate was estimated, so there is no estimate to assess.

`evidence_completeness` is `complete` when every alternative setting gave a rate and `incomplete` otherwise. The alternative settings include the entry `closure`, the rate after iteration, which is missing when the iteration did not converge. With incomplete evidence a provisional log can still be written. The missing settings are then listed in `report.json`, in `report.md`, in `warnings` of `provisional_manifest.json` and on the preview. Nothing is established about a missing estimate.

When the estimate is withheld, `run --export` writes no provisional log from it and states the reason:

```text
No provisional log was written from this acquisition estimate: <reason>. The estimate and failed checks remain in the report. A documented reference rate may be supplied.
```

Without a reference rate, the printed line adds `An uncorrected inspection view is available.` once the preview has been drawn. That preview shows the original or declared alignment without a rate correction. It draws the estimate in a separate panel whose title starts with `rejected candidate for inspection`. The estimate is not drawn as an applied correction, and `applied_mapping` records that no rate is in use. Without `--export` the report records the same status, and the preview follows it.

What you can do next:

- Read `reasons` and the numbers of the triggered condition, then the checks that are not met. They say why the estimate does not have a clear direction or a settled value.
- Declare a starting alignment that has an independent basis, or choose a matching signal with clearer contrast, as described above. Each is a new run with its own report.
- Supply a reference rate with `--reference-rate-ppm`. It is written as its own rate-only provisional log, `provisional_rate_source` is `user_reference_rate`, and `reference_rate.reason` is `acquisition_estimate_not_eligible`. The withheld estimate stays in the report and in the preview as the rejected candidate.

The eligibility depends on the tolerance of the run. A run with another tolerance is another assessment with its own record. Choose the tolerance for the intended use, not to obtain a log.

After a failed read-back the same eligibility applies to the estimate that the checked candidate was written from. An eligible candidate is published unchanged as the provisional log. A candidate that is not eligible stays in the run directory for diagnosis only, and a supplied reference rate is written as its own rate-only log. Changed bytes remain an execution error with exit code 1.

## A tolerance in micrometres, and distances in the report

`analysis.tolerance_s` has no default. You can give the tolerance as a distance along the scan instead: `"tolerance_um": 75` in `analysis`, or `--tolerance-um 75` on the command line. The program converts it to seconds once, before the diagnosis, with the largest scan speed of the valid firing intervals of the selected acquisition. The estimator, the checks and every read-back then use that value in seconds. `tolerance_resolution` in `report.json` records the request, the speed, its source and `resolved_tolerance_s`. The `config.json` of the run holds the resolved `tolerance_s` with the same record, so a rerun of that file and `iolite-export` do not convert again.

```text
python -B -X utf8 -m lagcheck run config.json --output ../trial_distance --tolerance-um 75
```

If the run stops with `A tolerance in micrometres cannot be converted to seconds for this record`, the message names the cause:

- The laser log gives no scan speed, or some firing intervals have no usable speed. A speed is not usable when its cell is blank with no value above it, when its magnitude is not a positive finite number, when it changes during the firing, or when the recorded positions and times of the firing interval contradict it by more than a factor of two. A speed that the firing start takes from a record above it is also not usable when its firing interval has no position to check it. The message gives the number of firing intervals for each reason. Give the tolerance in seconds, or declare the scan speed with `analysis.scan_speed_um_s`. A declared speed is used as given, and a difference from the speed in the log is recorded.
- The record contains spots: firing intervals without stage movement. A time difference has no equivalent distance along a scan there. Give the tolerance in seconds.
- The log is a New Wave log. Its speed field can hold a positioning speed, so no firing speed is taken from it and a declared speed is not applied to it. Give the tolerance in seconds.

Giving both `tolerance_s` and `tolerance_um` in the configuration is an error, and so is giving both `--tolerance-s` and `--tolerance-um`. One value on the command line replaces the tolerance of the configuration, whatever its unit, and the replaced value is kept in `tolerance_resolution`. Values that are not positive finite numbers are rejected.

Only the speed in force at each firing start is read: the speed cell of the On record, or the last value above it when the cell is blank. A speed written in an Off record that positions the stage is not a firing speed. The positions of a firing interval are used for that interval only. A speed that they contradict by more than a factor of two is set aside, and one firing interval without stage movement makes the record a spot record. A coordinate that is missing in one firing interval does not switch these checks off for the others. A firing interval without positions cannot be checked. A speed written in its On record is then used. A speed taken from a record above it is then not used, because that record can hold a positioning speed. The factor of two is a coarse screen: a speed that passes it is not shown to be correct. `distance_scale.speed` in `report.json` counts the firing intervals by where their speed was written (`log_speed_origins`), by the reason a speed was not used (`log_speeds_not_used`) and without positions (`intervals_without_positions`). With several scan speeds the largest one converts the request, so the tolerance in seconds is the strictest, and the report gives the range.

`distance_scale` in `report.json` and the section `Tolerance and equivalent distances` of `report.md` give time differences as distances along the scan when a scan speed is known: the tolerance, the largest held-out prediction errors, the range of the alternative settings and the residual rate over the firing span. A rate is multiplied by the firing span and by 10⁻⁶ before it is converted. Coverage fractions, correlation values, counts and results that are true or false have no equivalent distance. `Equivalent distances: not computable` means that no scan speed is known for the reason given; a run in seconds is complete without it. A distance is a scale along a line scanned at that speed. It does not establish the position error across a line change or a gap, and a logged speed is a nominal value.

Spot widths are given only when the width of the spot along the scan is known. Declare it with `analysis.spot_width_along_scan_um`. From the laser log it is taken only for a spot recorded as a circle, or as a square with its sides along the stage axes that is scanned along an axis. The recorded size must be a positive finite number and the angle of a square a finite number, and size, type and angle must not change while the laser fires. A size such as `15 x 15` without a recorded spot type is not interpreted. A spot width is a length scale. It is not an image resolution.

## Check what the new trial actually produced

Read the new report and preview metadata. Distinguish the input diagnostic, processed-input diagnostic when applicable, read-back, and export status. Inspect the stated failed checks. A provisional log stays provisional even if its preview looks better. Default plots are inspection views and may be a line-order/time display, a reduced grid, or a raw signal timeline; none establishes independent spatial accuracy or compatibility with a receiving importer.

If these explicit options do not yield a rate, retain the raw/declared-alignment view and report the limitation. The automatic retry preserves its stated input-placement assumption and provisional status. The program does not choose another acquisition or establish an unknown physical time origin.
