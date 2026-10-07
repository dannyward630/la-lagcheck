# Troubleshooting when no rate is available

Companion guide for LA-LagCheck 0.5.0. Run the commands from the package directory with its requirements installed. Keep configuration copies beside the original configuration so relative input paths retain their meaning.

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

## Explicitly choose a different matching signal

If the configured signal has little firing-pattern contrast, choose a recorded channel that has interpretable contrast for the same acquisition. Copy the configuration to `config_single_channel.json` and change the existing CSV `ms.signal_columns` list. For example, if the input contains an exact `Fe56` column:

```json
"signal_columns": ["Fe56"]
```

For a NetCDF mapping, change its existing `signal_variables` list instead, using an actual variable name from that file. Keep its time mapping and acquisition selection unchanged. A list of multiple channels is summed for diagnosis; a list of one channel makes that channel the matching signal.

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
