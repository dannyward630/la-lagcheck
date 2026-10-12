# Troubleshooting when no rate is available

Companion guide for LA-LagCheck 0.7.0. Run the commands from the package directory with its requirements installed. Keep configuration copies beside the original configuration so relative input paths retain their meaning.

Open `report.json` and `preview_pair.json` first. Read `reason`, `alignment`, `rate_ppm`, `warnings`, and the mapping source. A raw or unregistered view gives you something to inspect; it does not mean a correction was estimated. If `INCOMPLETE.json` exists, the run has an input or execution error. Its raw inspection plot does not make a corrected log usable.

Use a new configuration copy and a new output directory for each trial. Keep the original input files and previous results. The program first retains its strict diagnostic. With provisional output enabled, an automatic input-placement retry may supply a provisional rate under the conditions below. The later configuration edits are explicit user choices.

## Understand the automatic input-placement retry

If you have not declared a starting alignment, the program can repeat diagnosis when automatic alignment is unidentifiable or its rate fails the agreement or held-out prediction checks. It uses the parent input placement (`alignment_offset_s = 0` for an unsegmented input) with the same matching signal and other settings. Segment offsets preserve that placement by accounting for the shift from the parent origin; each segment is not reset to the start of the mass-spectrometer record. It does not choose a different channel. The retry is used only if the candidate rate is finite, the time-scale factor is positive, and agreement and prediction explicitly pass at both the initial-fit and iteration stages. Missing checks do not pass. Iteration convergence is reported but is not an additional requirement for provisional adoption. Its correction is provisional, and the strict result remains in the report.

Read the top-level `status`, `reason`, `rate_ppm`, `alignment`, `checks` and `input_diagnostic_status` as the unchanged strict result. Read `alignment_fallback` separately: it records `attempted`, `adopted`, `status`, `triggered_reason`, `alignment_offset_s`, `required_checks`, `failed_adoption_checks`, `diagnostic` and `note`. The candidate diagnostic is also saved in `alignment_fallback_report.json` and its block files. Adoption is marked `provisional_adopted` with rate source `input_position_fallback`; the run returns exit code 2. Failed agreement or prediction checks remain visible.

To disable this automatic retry while retaining other provisional outputs, run:

```text
python -B -X utf8 -m lagcheck run config.json --output ../trial_without_fallback --export --no-alignment-fallback
```

The same switch is available in `wizard`. `--no-provisional` also disables the retry; `--no-preview` only disables plots. A command without `--export` may still run the retry and show its provisional mapping, but writes no corrected log. An adopted retry takes precedence over the original estimate and an explicitly supplied reference. If the retry is rejected, an original rate can still be shown or exported under the existing provisional rules, with both the original failed checks and the rejected retry checks recorded.

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

## Supply a measured reference rate when the acquisition has no estimate

If a rate measured for the same instrument pair is appropriate and you elect to use it, supply it explicitly:

```text
python -B -X utf8 -m lagcheck run config.json --output ../trial_reference --export --reference-rate-ppm 170.9
```

`170.9` is the study's instrument-pair example, not a universal default. Replace it with your justified reference value and retain its source. An adopted input-placement retry takes precedence over both the original strict estimate and a supplied reference. The reference is used only when no acquisition estimate is available after this selection. It does not replace an existing estimate, supply a missing time origin, repair invalid input, or turn an unsupported diagnostic into a supported one.

When input and export requirements are satisfied, this route can produce `provisional_laser_log.csv` and `provisional_manifest.json`, marked not verified with the reference rate source. Without a valid starting alignment, fixed-mapping read-back can remain `not_run`; the display says that its clock origins are unregistered. Do not treat the plot as confirmation of absolute alignment.

## Check what the new trial actually produced

Read the new report and preview metadata. Distinguish the input diagnostic, processed-input diagnostic when applicable, read-back, and export status. Inspect the stated failed checks. A provisional log stays provisional even if its preview looks better. Default plots are inspection views and may be a line-order/time display, a reduced grid, or a raw signal timeline; none establishes independent spatial accuracy or compatibility with a receiving importer.

If these explicit options do not yield a rate, retain the raw/declared-alignment view and report the limitation. The automatic retry preserves its stated input-placement assumption and provisional status. The program does not choose another acquisition or establish an unknown physical time origin.
