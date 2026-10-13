# Quick start

The program needs Python 3.10 or newer. This release was validated with Python 3.12.5 and the pinned versions of `requirements-reference.txt`, which have no prebuilt packages for Python 3.13 or later (see `README.md`). Run from the package directory with dependencies installed as shown in `README.md`. Use a new output directory for each run. Keep `-X utf8` in the commands: the program prints µ and superscript characters, and without it a redirected output on Windows can fail at the last printed line.

## Chromium-grating reproduction

The [grating README](examples/chromium_grating/README.md) provides two portable configurations, raw files, commands and expected results for manuscript Section 3.4. Its declared start is specific to that acquisition.

## Synthetic geometry and iolite ESL export

This example has 240 scan lines, zero rate drift and a known constant offset. The first stored MS sample is defined by the generator.

```text
python -B -X utf8 examples/make_iolite_demo.py --output ../demo_inputs
python -B -X utf8 -m lagcheck run ../demo_inputs/config.json --output ../demo_canonical --export --export-mode affine --ms-first-sample-time "2026-01-01 00:00:00.000000" --reference-note "Known synthetic first stored sample"
python -B -X utf8 -m lagcheck iolite-export ../demo_canonical --output ../demo_iolite
```

`iolite-export` checks formatting and timing without opening iolite. Retain the canonical output directory and its original inputs for this second step. Verify the exact exported file in the receiving application before image use.

Expected output: `run` prints `"status": "supported"`, `"rate_ppm": 0.0` and `"correction_exported": true`, and its exit code is 0. The generator's `truth.json` gives the same rate, 0.0, and a physical offset of 3 s. In `report.json`, `post_export_verification` has `passed` true and a residual rate of 0.0. `iolite-export` prints `"status": "internal_checks_passed_importer_validation_pending"`, its exit code is 0, and it writes `corrected_laser_log_iolite.csv`. Compare numbers with an absolute tolerance of 0.000001, the tolerance of the grating verifier.

## Synthetic rate drift

```text
python -B -X utf8 examples/make_demo.py --output ../demo_rate_inputs
python -B -X utf8 -m lagcheck run ../demo_rate_inputs/config.json --output ../demo_rate_canonical --export
```

The generator writes the known rate and offset in its truth JSON. This interval-only example lacks the stage geometry required for iolite ESL export.

Expected output: `run` prints `"status": "supported"`, `"rate_ppm": 100.69370075580798`, `"rate_kind": "iterated"` and `"correction_exported": true`, and its exit code is 0. The iteration converges at its first trial rate, so the rate after iteration equals the single-pass rate here, and `"single_pass_rate_se_ppm"` is 1.408494111994502. The generator's `truth.json` gives a rate of 100.0. The estimate differs from it by 0.69 microseconds per second, which is less than the standard error of the single-pass slope. The package tests use the same synthetic record and accept a difference of up to 3 microseconds per second from 100. To compare a rerun with the numbers listed here, use an absolute tolerance of 0.000001. The last digits can differ between systems.

## Configure a new acquisition

Configuration paths are relative to the configuration file's directory. Use the exact laser sample selector for the intended acquisition. For example:

```json
{
  "ms": {
    "path": "signal.csv",
    "format": "agilent_masshunter_csv",
    "signal_columns": ["C13", "Fe56"]
  },
  "laser": {
    "path": "laser.csv",
    "format": "laser_state_csv",
    "time_kind": "datetime",
    "sample_selector": "Exact sample name"
  },
  "analysis": {"tolerance_s": 0.2},
  "export": {"mode": "rate_only"}
}
```

`ms.signal_columns` must name columns of your signal file; `C13` and `Fe56` are the channels of this example. If a name is not found, the error message lists the columns of the file. For `laser_state_csv` and `new_wave_laser_log`, `laser.sample_selector` is compared with the `Comment` column of the laser log, or with the column named by `laser.sample_column`. A blank cell takes the last nonblank value above it. The selector must equal that text exactly. The comparison is case-sensitive, and spaces around the cell text are ignored. If the log holds one sample name only, the selector can be left out. If the selector is not found, the error message lists the sample names of that column. `tolerance_s` has no default and must be given; the example uses 0.2 s.

The tolerance can be given as a distance along the scan instead of seconds: write `"tolerance_um": 75` in place of `"tolerance_s"`, or add `--tolerance-um 75` to the command. Give one of the two, not both. The program converts the distance to seconds once, before the diagnosis, with the largest scan speed of the firing intervals. For a log that records 375 µm s⁻¹, 75 µm is 0.2 s, and the result is that of `"tolerance_s": 0.2`. `tolerance_resolution` in `report.json` records the request, the speed and its source. `--tolerance-s` and `--tolerance-um` replace the tolerance of the configuration for one run. A distance needs a scan speed for every firing interval: from the laser log, or declared with `"scan_speed_um_s"`. It is not accepted for spot records and New Wave logs, which need seconds. See [TROUBLESHOOTING.md](TROUBLESHOOTING.md) when a distance is refused.

With several `signal_columns` in a CSV configuration, the program compares the sum and each column by the peak of the automatic starting-alignment search and uses the highest (`ms.signal_selection`, default `"highest_alignment_peak"`; `"sum"` uses the sum without comparison). If a block peak lies at a search boundary imposed by the start or end of the signal record, it repeats the diagnosis without the scan lines within 30 s of that end (`analysis.skip_record_edge_lines`, default `true`; see `record_edge_lines` in the report). Both defaults were set after examining a set of existing records. The example tolerance is a run setting, not a general accuracy guarantee. See the other JSON examples for Thermo and PerkinElmer reader mappings. The default sensitivity analysis requires at least 180 firing intervals (15 blocks with at least 12 intervals each); more may be needed if the baseline block count exceeds 15. This count alone does not establish sufficient timing evidence. The count that matters is that of the firing intervals in the fits, after any scan lines at the ends of the record have been left out. With fewer, the report says so in `notes` and names the stability alternatives that could not be estimated.

Raw readers rebase each time axis to its first retained timestamp. A declared offset must refer to these rebased coordinates, not to the difference between calendar headers. Automatic alignment searches for a start and checks its uniqueness. An unsupported decision does not establish damaged physical data or require discarding it; inspect the reasons and establish alignment for reconstruction.

## Export interpretation

A rate-only export applies `t' = t0 + (1 + rate_ppm * 1e-6) * (t - t0)`. Determine its constant offset in the reconstruction importer. Both firing-interval endpoints and continuous-trajectory transitions must receive the same time-coordinate mapping.

An affine export additionally needs the first stored MS sample's time in the original laser-log coordinate and a provenance note (`--ms-first-sample-time`, `--reference-note`). An acquisition header is not necessarily that first-sample time. Do not apply the original rate again to either corrected log.

For `run` and `iolite-export`, exit code 0 means the requested internal checks passed; 2 means a diagnostic or export check did not pass; 1 means an input or execution error. Invalid command-line syntax also returns 2. For `screen-log`, 0 means that at least one sequence was screened and no interval was flagged; 2 means a flagged interval or that nothing could be screened, and `status` in the printed result distinguishes the two; `unevaluated_segments` counts segments with fewer than three line-change intervals. Inspect `report.json`, including `status`, `correction_exported`, `usable_file`, and the export verification. A completed canonical export has a `correction_manifest.json` identifying its `output_file` and `output_sha256`; do not treat files in an incomplete output directory as accepted corrections.

## What the report states

`Rate estimate` in `report.md` is `rate_ppm` in `report.json` and in the line that `run` prints at the end. It is the rate after iteration when the iteration converged; the printed line then has `"rate_kind": "iterated"`. When the iteration did not converge, it is the single-pass rate (`"rate_kind": "single_pass"`), and a note gives the last trial rate and the residual rate left at that trial.

The single-pass rate is always stated together with the standard error of its slope (`baseline.rate_ppm` and `baseline.slope_se_ppm` in `report.json`; `single_pass_rate_ppm` and `single_pass_rate_se_ppm` in the printed line). This is the ordinary least-squares standard error of the slope of block lag against time over the retained blocks. It is conditional on those blocks and on a single line. It belongs to the single-pass rate. The rate after iteration is stated without a standard error.

For each check that is not met, `report.md` gives the measured value and its limit where the report holds them. The numbers come from `report.json` and `blocks.csv`. The held-out error of a stability alternative is not stored, so for an alternative that fails on it the text names that condition without a number. When the diagnosis stops before the checks (`insufficient_evidence`), the section gives the reason with its numbers and, when blocks were measured, points to `blocks.csv`.

`notes` in `report.json` is repeated in `report.md` and in the printed line. It says when the record has too few firing intervals for a fit or for every stability alternative, when the iteration did not converge, when no provisional log was written from the estimate, and whether a reference rate given with `--reference-rate-ppm` was used. The printed line adds that an uncorrected inspection view is available once the preview has been drawn. `reference_rate` in `report.json` holds the same statement with a `used` flag and a `reason`. Without `--export`, `report.md` says how to request a corrected or provisional log.

`provisional_export_eligibility` in `report.json` says whether the rate estimated from the acquisition may be written to a provisional log. `status` is `not_blocked`, `withheld`, `not_assessable` or `no_acquisition_rate`. `evidence_completeness` is `complete`, or `incomplete` when an alternative setting gave no rate. The object holds the candidate rate and its source, the single-pass rate with the standard error of its slope and the number of blocks, the firing span, the tolerance, the numbers of the three conditions and the mapping in use (`applied_mapping`). `not_blocked` means only that none of the three conditions was triggered. `report.md` repeats the object in the section `Provisional export eligibility`. The section `Tolerance and equivalent distances` states the tolerance of the run and how it was given. When a scan speed is known, it gives the time differences of the diagnosis as distances along the scan (`distance_scale` in `report.json`).

## Automatic input-placement retry

When starting alignment is not declared, an unidentifiable automatic start can trigger a separate retry at the original input placement. An automatic rate that fails block agreement or held-out prediction (including their iteration-stage checks) also triggers this retry. It uses the same signal and settings, setting the input-placement offset to zero, with each fragment's original origin shift retained. The candidate is adopted only with a finite rate, positive time scale, and explicitly passing agreement and prediction at both initial-fit and iteration stages. Convergence is reported but is not an additional adoption gate. Top-level strict diagnostic fields stay unchanged; the candidate is recorded under `alignment_fallback` and in `alignment_fallback_report.json` and its block files. Adoption is always provisional, even if read-back passes, and exits 2. Adopted candidates take precedence over the original estimate and a reference rate. An adopted candidate is written to a log only when it is eligible for provisional export; `provisional_export_eligibility` is assessed for the adopted candidate, not for the strict diagnostic. A rejected retry retains the original candidate or reference under the existing rules and records its rejection. `--no-alignment-fallback` disables retry; `--no-provisional` disables retry and provisional logs; `--no-preview` only disables drawing. Without `--export`, retry can still supply a preview but no log. See [TROUBLESHOOTING.md](TROUBLESHOOTING.md) for the input-placement assumption and explicit next steps.

## Provisional output when the checks do not pass

With `--export`, a correction is labelled verified only when the original input requires no zero-duration exclusion or automatic segmentation, and its input diagnostic and fixed-mapping read-back both pass. When the timing checks fail at either stage and a rate is available and eligible for provisional export, the run writes `provisional_laser_log.csv` and, last, `provisional_manifest.json` (`provisional: true`, `verified: false`, `failure_stage` = `input_diagnostic`, `alignment_fallback` or `readback`, the failed checks, the rate and its source, `readback_status`, `provisional_export_eligibility`, `tolerance_resolution` and `warnings`). After a failed read-back, the provisional file keeps the already checked candidate bytes when the estimate is eligible for provisional export; the correction is not applied again. When the estimate is not eligible, the candidate stays in the run directory for diagnosis only, and a reference rate, if one was given, is written as its own rate-only provisional log. If those bytes changed after they were checked, the run stops as an execution error (exit code 1) and writes no provisional file. After a failed input diagnostic, the rate is the one estimated from the acquisition; when none can be estimated, or when the estimate is not eligible, give a reference rate measured on the same instrument pair, for example `--reference-rate-ppm 170.9`. Without an eligible estimate or a reference rate, no provisional log is written. An affine provisional log is written when an affine export is requested and the fitted constant of a converged iteration is available; otherwise the provisional log is rate-only and the report says why. The written bytes, rate and constant are checked against the stated mapping before publication. When a starting alignment is available, a provisional log that was written as its own file is also re-diagnosed with the mapping held fixed (`provisional_readback_report.json`); otherwise `readback_status` is `not_run` and no new alignment is guessed. A provisional file that keeps the candidate bytes has the read-back of that candidate (`post_export_report.json`) and is not diagnosed again. Export configuration errors (unknown keys, an unknown mode, an unreadable time reference, an invalid reference rate) are input errors with exit code 1, and no corrected log is published. A provisional file, or its iolite copy, is refused as input to a further `--export`, so a rate cannot be applied twice. `iolite-export` formats a provisional run in the same way and writes `provisional_laser_log_iolite.csv` with its own `provisional_manifest.json`. A provisional file is not a verified correction and does not verify the absolute alignment: inspect the reconstruction and state the failed stage and checks with any result. Use `--no-provisional` to write only verified corrections.

An estimate is not eligible in three cases: its direction is uncertain; the alternative settings give rates of both signs; the iteration did not converge and the rate still changes. In each case the effect over the firing span must exceed the tolerance. The run then writes no provisional log from the estimate and says why: `No provisional log was written from this acquisition estimate: <reason>. The estimate and failed checks remain in the report. A documented reference rate may be supplied.` The exit code stays 2, and the estimate and the checks stay in the report. A reference rate is written as a rate-only log that takes no fitted constant from the estimate. [TROUBLESHOOTING.md](TROUBLESHOOTING.md) gives the three conditions, what a provisional log can be used for and which numbers to read before using it.

`run` and `wizard` draw `preview_pair.png` and `preview_pair.json` by default, including when timing checks fail. The JSON records the actual displayed rate and its source, input and read-back checks, publication status, coordinate provenance, warnings and any grid reduction. Use `--no-preview` to disable this explicitly, or `--preview Fe56` to choose a channel. `--no-sample-label` leaves the sample label of the laser log out of the title of `preview_pair.png`. When the log has no sample label, the title shows the name of the laser file instead, and the switch leaves that out too. Without the switch the title is unchanged. `preview` and `wizard` accept the same switch. The additional New Wave views keep the sample labels of their panels. A later `python -B -X utf8 -m lagcheck preview <result-directory> --name another_view` draws another view without overwriting outputs. Its exit 0 means the view was drawn, not that the timing or spatial result passed.

A figure based on a published rate-only log uses the fixed diagnostic alignment and says so; it is not a tested importer offset. Without a known alignment, zero is only a display shift and is labelled unregistered. A reference rate is used only when the acquisition has no estimate or its estimate is not eligible for provisional export. No laboratory rate is silently inserted. An estimate that is withheld from the provisional log is not drawn as an applied correction. It then appears in a separate panel whose title starts with `rejected candidate for inspection`, beside the original or declared alignment and, when a reference rate was supplied, beside that rate. For a result that is not supported, the footer of the figure states the tolerance, the eligibility and the alternative settings that gave no rate. No coordinate data means a line-order/time view, not fabricated spatial geometry. If no samples land in firing intervals, the raw signal timeline shows what is available and lists the reason.

The default zero-duration policy excludes zero-duration pairs from analysis while retaining and reporting their source records. To require the old strict reader behavior, set `laser.zero_duration_policy` to `reject`. Original inputs are never edited. Corrected CSVs retain those event rows, and their read-back configuration carries the exclusion policy; receiving-importer handling is not established by this program. For the generic laser-state adapter, selected-sample fragments separated by other samples or sequence changes get separate subdirectories, reports and previews; inspect `segments.json` for every fragment, empty metadata fragment, source hash and origin shift. No joint correction is exported. A malformed fragment is reported as an error even when another valid fragment can be inspected.

The report distinguishes the original input check from the processed-interval diagnostic. Logs made after exclusion or automatic segmentation are provisional even if the processed diagnostic and read-back pass. Input and byte-integrity errors remain exit 1; an unchanged readable MS file may still receive a raw diagnostic view. A failed requested rendering also returns 1 and records `preview_error.json`, even if that fallback view succeeds.

For explicitly mapped numeric CSV inputs with missing/nonfinite signal cells, an input-error run can draw the finite raw signal and mark gaps, with every affected original data-record number in `input_processing.json`. This is inspection only: exit 1 and `INCOMPLETE.json` remain, no cleaned series is analysed, and no correction is published.


For New Wave event logs, use the explicit adapter and selection shown in [NEW_WAVE.md](NEW_WAVE.md). Repeated On records and changing line/spot sequence numbers are supported under those semantics; conditional short-record results retain warnings and do not become verified corrections.

The explicit New Wave adapter has different documented sequence semantics: line/spot numbering remains within the selected acquisition. An exact sample selector may omit interleaved standards without bridging any firing interval; see the retained source-record gaps in `input_processing.json`.
