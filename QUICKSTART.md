# Quick start

Use Python 3.10 or newer; the pinned reference environment was tested with Python 3.12.7. Run from the package directory with dependencies installed as shown in `README.md`. Use a new output directory for each run.

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

## Synthetic rate drift

```text
python -B -X utf8 examples/make_demo.py --output ../demo_rate_inputs
python -B -X utf8 -m lagcheck run ../demo_rate_inputs/config.json --output ../demo_rate_canonical --export
```

The generator writes the known rate and offset in its truth JSON. This interval-only example lacks the stage geometry required for iolite ESL export.

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

The example tolerance is a run setting, not a general accuracy guarantee. See the other JSON examples for Thermo and PerkinElmer reader mappings. The default sensitivity analysis requires at least 180 firing intervals (15 blocks with at least 12 intervals each); more may be needed if the baseline block count exceeds 15. This count alone does not establish sufficient timing evidence.

Raw readers rebase each time axis to its first retained timestamp. A declared offset must refer to these rebased coordinates, not to the difference between calendar headers. Automatic alignment searches for a start and checks its uniqueness. An unsupported decision does not establish damaged physical data or require discarding it; inspect the reasons and establish alignment for reconstruction.

## Export interpretation

A rate-only export applies `t' = t0 + (1 + rate_ppm * 1e-6) * (t - t0)`. Determine its constant offset in the reconstruction importer. Both firing-interval endpoints and continuous-trajectory transitions must receive the same time-coordinate mapping.

An affine export additionally needs the first stored MS sample's time in the original laser-log coordinate and a provenance note (`--ms-first-sample-time`, `--reference-note`). An acquisition header is not necessarily that first-sample time. Do not apply the original rate again to either corrected log.

For `run` and `iolite-export`, exit code 0 means the requested internal checks passed; 2 means a diagnostic or export check did not pass; 1 means an input or execution error. Invalid command-line syntax also returns 2. Inspect `report.json`, including `status`, `correction_exported`, `usable_file`, and the export verification. A completed canonical export has a `correction_manifest.json` identifying its `output_file` and `output_sha256`; do not treat files in an incomplete output directory as accepted corrections.

## Automatic input-placement retry

When starting alignment is not declared, an unidentifiable automatic start can trigger a separate retry at the original input placement. An automatic rate that fails block agreement or held-out prediction (including their iteration-stage checks) also triggers this retry. It uses the same signal and settings, setting the input-placement offset to zero, with each fragment's original origin shift retained. The candidate is adopted only with a finite rate, positive time scale, and explicitly passing agreement and prediction at both initial-fit and iteration stages. Convergence is reported but is not an additional adoption gate. Top-level strict diagnostic fields stay unchanged; the candidate is recorded under `alignment_fallback` and in `alignment_fallback_report.json` and its block files. Adoption is always provisional, even if read-back passes, and exits 2. Adopted candidates take precedence over the original estimate and a reference rate. A rejected retry retains the original candidate or reference under the existing rules and records its rejection. `--no-alignment-fallback` disables retry; `--no-provisional` disables retry and provisional logs; `--no-preview` only disables drawing. Without `--export`, retry can still supply a preview but no log. See [TROUBLESHOOTING.md](TROUBLESHOOTING.md) for the input-placement assumption and explicit next steps.

## Provisional output when the checks do not pass

With `--export`, a correction is labelled verified only when the original input requires no zero-duration exclusion or automatic segmentation, and its input diagnostic and fixed-mapping read-back both pass. When the timing checks fail at either stage and a rate is available, the run writes `provisional_laser_log.csv` and, last, `provisional_manifest.json` (`provisional: true`, `verified: false`, `failure_stage` = `input_diagnostic`, `alignment_fallback` or `readback`, the failed checks, the rate and its source, and `readback_status`). After a failed read-back, the provisional file keeps the already checked candidate bytes; the correction is not applied again. If those bytes changed after they were checked, the run stops as an execution error (exit code 1) and writes no provisional file. After a failed input diagnostic, the rate is the one estimated from the acquisition; when none can be estimated, give a reference rate measured on the same instrument pair, for example `--reference-rate-ppm 170.9`. Without either rate, no log is written. An affine provisional log is written when an affine export is requested and the fitted constant of a converged iteration is available; otherwise the provisional log is rate-only and the report says why. The written bytes, rate and constant are checked against the stated mapping before publication. When a starting alignment is available, the written file is also re-diagnosed with the mapping held fixed (`provisional_readback_report.json`); otherwise `readback_status` is `not_run` and no new alignment is guessed. Export configuration errors (unknown keys, an unknown mode, an unreadable time reference, an invalid reference rate) are input errors with exit code 1, and no corrected log is published. A provisional file, or its iolite copy, is refused as input to a further `--export`, so a rate cannot be applied twice. `iolite-export` formats a provisional run in the same way and writes `provisional_laser_log_iolite.csv` with its own `provisional_manifest.json`. A provisional file is not a verified correction and does not verify the absolute alignment: inspect the reconstruction and state the failed stage and checks with any result. Use `--no-provisional` to write only verified corrections.

`run` and `wizard` draw `preview_pair.png` and `preview_pair.json` by default, including when timing checks fail. The JSON records the actual displayed rate and its source, input and read-back checks, publication status, coordinate provenance, warnings and any grid reduction. Use `--no-preview` to disable this explicitly, or `--preview Fe56` to choose a channel. A later `python -B -X utf8 -m lagcheck preview <result-directory> --name another_view` draws another view without overwriting outputs. Its exit 0 means the view was drawn, not that the timing or spatial result passed.

A figure based on a published rate-only log uses the fixed diagnostic alignment and says so; it is not a tested importer offset. Without a known alignment, zero is only a display shift and is labelled unregistered. A reference rate is used only when the acquisition has no estimate. No laboratory rate is silently inserted. No coordinate data means a line-order/time view, not fabricated spatial geometry. If no samples land in firing intervals, the raw signal timeline shows what is available and lists the reason.

The default zero-duration policy excludes zero-duration pairs from analysis while retaining and reporting their source records. To require the old strict reader behavior, set `laser.zero_duration_policy` to `reject`. Original inputs are never edited. Corrected CSVs retain those event rows, and their read-back configuration carries the exclusion policy; receiving-importer handling is not established by this program. For the generic laser-state adapter, selected-sample fragments separated by other samples or sequence changes get separate subdirectories, reports and previews; inspect `segments.json` for every fragment, empty metadata fragment, source hash and origin shift. No joint correction is exported. A malformed fragment is reported as an error even when another valid fragment can be inspected.

The report distinguishes the original input check from the processed-interval diagnostic. Logs made after exclusion or automatic segmentation are provisional even if the processed diagnostic and read-back pass. Input and byte-integrity errors remain exit 1; an unchanged readable MS file may still receive a raw diagnostic view. A failed requested rendering also returns 1 and records `preview_error.json`, even if that fallback view succeeds.

For explicitly mapped numeric CSV inputs with missing/nonfinite signal cells, an input-error run can draw the finite raw signal and mark gaps, with every affected original data-record number in `input_processing.json`. This is inspection only: exit 1 and `INCOMPLETE.json` remain, no cleaned series is analysed, and no correction is published.


For New Wave event logs, use the explicit adapter and selection shown in [NEW_WAVE.md](NEW_WAVE.md). Repeated On records and changing line/spot sequence numbers are supported under those semantics; conditional short-record results retain warnings and do not become verified corrections.

The explicit New Wave adapter has different documented sequence semantics: line/spot numbering remains within the selected acquisition. An exact sample selector may omit interleaved standards without bridging any firing interval; see the retained source-record gaps in `input_processing.json`.
