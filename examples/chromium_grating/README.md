# Chromium-grating example

This example contains the MassHunter signal CSV with its local instrument-path header removed, and the original laser log from a chromium grating scanned in 621 lines over 8.8 h. Si29 is the matching signal. It reproduces the raw-file timing diagnostic and saved-log readback in manuscript section 3.4 and SI Sections S9 and S12. The spatial comparison in the manuscript uses the automatic iterated rate of 173.11 µs s⁻¹. This example verifies the timing analysis; it does not reproduce the separate grating geometry measurements.

From the package root, install the dependencies if needed (`requirements-reference.txt` holds the validated versions, as in the main `README.md`; `requirements.txt` gives ranges), then run:

```text
python -B -X utf8 -m lagcheck run examples/chromium_grating/auto.json --output ../grating_results/automatic
python -B -X utf8 -m lagcheck run examples/chromium_grating/first_si_event.json --output ../grating_results/declared --export --export-mode rate_only
python -B -X utf8 examples/chromium_grating/verify_results.py --automatic ../grating_results/automatic --declared ../grating_results/declared
python -B -X utf8 examples/chromium_grating/verify_event_alignment.py
```

Each run requires a new output directory. Results are written outside the package. Configuration input paths are resolved relative to each configuration file; no original computer path is needed.

| Starting alignment | Iterated rate (µs s⁻¹) | Six checks | Decision | Export |
| --- | ---: | --- | --- | --- |
| Automatic | 173.11258783602625 | All pass; no line-change interval flagged | supported | Not requested |
| Declared first Si29 event | 173.05804299808756 | All pass; no line-change interval flagged | supported | Rate-only CSV, saved-file readback passes |

Starting alignment is assessed before the six checks: coverage, agreement among usable blocks, distance from search boundaries, robustness to alternative settings, held-out prediction and iteration. `expected_results.json` maps these to the individual Boolean fields in `report.json` and stores the expected numerical fields. The verifier checks these values within an absolute numerical tolerance of 0.000001, packaged-input hashes, and agreement of candidate, diagnosed and final CSV bytes. The last digits of these numbers can differ from one environment to another, so the comparison uses a tolerance and not an exact match. When this version was validated, the rate of the declared branch was 173.05804299808747, and the recorded value is 173.05804299808756. Its exit codes are 0 for a match, 2 for a mismatch and 1 for an input or tool error. The first line starts 22.68 s after the start of the record, so no block peak lies at a record-limited boundary and `record_edge_lines.applied` is false in both runs; the rates are identical to those of 0.6.0 within rounding. Full run details remain in each output directory.

The native MS parser subtracts the first recorded sample time, 0.0814 s. The laser parser uses the first metadata timestamp, placing first firing On at 23.637 s. The declared event uses the geometric mean of median Si29 baseline at raw MS time below 10 s and plateau at 25–40 s, then the first upward crossing before 45 s with linear interpolation. The crossing is 21.43557327577439 s on the original MS time axis, or 21.354173275774393 s after rebasing. This gives the fixed offset −2.2828267242256075 s. `coordinate_reference.json` records the calculation; `verify_event_alignment.py` recalculates it without changing the settings. A declared starting alignment is supplied evidence, not an independent check of the physical timing reference.

The declared branch produces `corrected_laser_log.csv`, `correction_manifest.json`, and `post_export_report.json`. The residual lag rate is 0.12545095339121373 µs s⁻¹, rounded to 0.13 in the paper. Saved-file readback retains the timing mapping and performs no automatic realignment. Candidate and final CSV bytes are identical. A constant import offset is still required for this rate-only log. These checks establish internal timing consistency of the saved file. This example does not test those bytes in a receiver GUI, reconstruct an image, or establish independent spatial accuracy.

Both configurations give the tolerance in seconds (0.2 s) and declare the scan speed of 375 µm s⁻¹. `distance_scale` in each `report.json` therefore gives that tolerance as 75 µm along the scan, and the other time differences of the diagnosis in micrometres. At every firing start of the laser log the speed in force is 375 µm s⁻¹. It is set in a record above the firing start and agrees with the recorded positions of each firing interval. The 20000 µm s⁻¹ of the positioning records is not read as a firing speed. The log gives the spot size as `15 x 15` without a spot type, so no spot widths are reported.

## The two affine logs

SI Section S9 also reports the read-back of two affine logs. They are written with `--export-mode affine`:

```text
python -B -X utf8 -m lagcheck run examples/chromium_grating/auto.json --output ../grating_results/automatic_affine --export --export-mode affine --ms-first-sample-time "2026-09-25 13:29:36.0814" --reference-note "MassHunter header acquisition time plus the first stored sample time"
python -B -X utf8 -m lagcheck run examples/chromium_grating/first_si_event.json --output ../grating_results/declared_affine --export --export-mode affine --ms-first-sample-time "2026-09-25 13:29:36.0814" --reference-note "MassHunter header acquisition time plus the first stored sample time"
```

`--ms-first-sample-time` declares the time of the first stored MS sample in the timestamp coordinate of the written log. The value used here is the acquisition time in the header of the signal file (`Acquired : 25-Sep-26 13:29:36`) plus the time of the first stored sample (0.0814 s). The program cannot verify this reference, and this example does not test it in a receiving application. The read-back subtracts the reference again before the comparison, so the read-back numbers do not depend on it. A run of the automatic branch with `2026-09-25 13:25:23.000000` gave the same read-back numbers.

| Log | Exit code | Status | Read-back residual rate (µs s⁻¹) | Coverage fraction | Largest held-out error (s) | Written log |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| Automatic, affine | 0 | supported | 0.10442991032827022 | 0.999540784259613 | 0.005658407779361066 | `corrected_laser_log.csv` |
| Declared first Si29 event, affine | 2 | export_qc_failed | 0.20926164006378908 | 0.9995407842345649 | 0.0049945499769129984 | `provisional_laser_log.csv` |

The numbers are in `post_export_verification` of each `report.json` (`residual_rate_ppm`, `coverage_fraction`, `heldout_max_s`). They round to the 0.104430 and 0.209262 µs s⁻¹ of SI Section S9. In the declared branch the residual rate exceeds the limit of 0.2 µs s⁻¹, the value of `closure_tolerance_ppm` in this run (`residual_rate_below_threshold` is not met). The input diagnostic of that run is supported, its read-back fails, and the log is published as provisional (`failure_stage` is `readback` in `provisional_manifest.json`). The estimate that the candidate was written from is eligible for provisional export (`provisional_export_eligibility.status` is `not_blocked` in `report.json`), so the checked candidate bytes are published unchanged. `verify_results.py` checks the runs of the first two commands of this file only. Compare these two runs with an absolute tolerance of 0.000001.

## Records that are not in time order

The preview footers and `coordinate_reference.json` contain the line `1241 decreasing metadata timestamps retained in original row order; no repair applied`. It means that 1241 records of this laser log carry a timestamp earlier than that of the record before them: the log is not written in time order. The program does not sort the records. It reads the firing intervals from the `On` and `Off` records in file order and stops with an error if those intervals are out of order or overlap. The line reports the count and is not an error. A corrected log keeps the row order of the original log. `iolite-export` writes the records in time order and records the correspondence in `row_mapping.csv` (see [IMPORTER_CHECK_INPUTS.md](../../IMPORTER_CHECK_INPUTS.md)).
