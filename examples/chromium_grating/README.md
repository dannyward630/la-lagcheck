# Chromium-grating example

This example contains the MassHunter signal CSV with its local instrument-path header removed, and the original laser log from a chromium grating scanned in 621 lines over 8.8 h. Si29 is the matching signal. It reproduces the raw-file timing diagnostic and saved-log readback in manuscript section 3.4 and SI Sections S9 and S12. The spatial comparison in the manuscript uses the automatic iterated rate of 173.11 µs s⁻¹. This example verifies the timing analysis; it does not reproduce the separate grating geometry measurements.

From the package root, install the dependencies in `requirements.txt` if needed, then run:

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

Starting alignment is assessed before the six checks: coverage, agreement among usable blocks, distance from search boundaries, robustness to alternative settings, held-out prediction and iteration. `expected_results.json` maps these to the individual Boolean fields in `report.json` and stores the expected numerical fields. The verifier checks these values within an absolute numerical tolerance of 0.000001, packaged-input hashes, and agreement of candidate, diagnosed and final CSV bytes. Its exit codes are 0 for a match, 2 for a mismatch and 1 for an input or tool error. The first line starts 22.68 s after the start of the record, so no block peak lies at a record-limited boundary and `record_edge_lines.applied` is false in both runs; the rates are identical to those of 0.6.0 within rounding. Full run details remain in each output directory.

The native MS parser subtracts the first recorded sample time, 0.0814 s. The laser parser uses the first metadata timestamp, placing first firing On at 23.637 s. The declared event uses the geometric mean of median Si29 baseline at raw MS time below 10 s and plateau at 25–40 s, then the first upward crossing before 45 s with linear interpolation. The crossing is 21.43557327577439 s on the original MS time axis, or 21.354173275774393 s after rebasing. This gives the fixed offset −2.2828267242256075 s. `coordinate_reference.json` records the calculation; `verify_event_alignment.py` recalculates it without changing the settings. A declared starting alignment is supplied evidence, not an independent check of the physical timing reference.

The declared branch produces `corrected_laser_log.csv`, `correction_manifest.json`, and `post_export_report.json`. The residual lag rate is 0.12545095339121373 µs s⁻¹, rounded to 0.13 in the paper. Saved-file readback retains the timing mapping and performs no automatic realignment. Candidate and final CSV bytes are identical. A constant import offset is still required for this rate-only log. These checks establish internal timing consistency of the saved file. This example does not test those bytes in a receiver GUI, reconstruct an image, or establish independent spatial accuracy.
