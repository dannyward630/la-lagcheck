# Receiving-application checks

`iolite-export` formats a verified canonical laser log as an ESL CSV and checks millisecond timestamps, stable event ordering, firing-event pairing, unchanged non-time fields and the declared fixed time mapping. It does not open iolite or HDIP, save a reconstruction workspace, or validate an image in either application.

`iolite-export` writes the records in time order. It rounds each timestamp to the millisecond and sorts the records by the rounded value. Records with the same timestamp keep their original order. A laser log whose records are not in time order therefore has a different row order after this step. The canonical log (`corrected_laser_log.csv` or `provisional_laser_log.csv`) keeps the row order of the original log. `format_stage/row_mapping.csv` gives the correspondence. It has one row per output record, with the number of the record in the output (`output_record`), its number in the canonical log (`source_record`), both timestamps, the rounding difference and the role of the record (`firing_start`, `firing_end` or `metadata`). When sorting would change the sample or sequence that a record inherits from the rows above it, the step writes the inherited value into the blank cell and lists it in `format_stage/compatibility_receipt.json` (`materialized_membership_cells`). The step stops with an error if sorting would change the pairs of firing endpoints or the scan velocity at a firing endpoint.

Keep the original inputs, configuration, canonical export, receipts and SHA256 hashes together. The original input paths must remain accessible when the ESL formatting step rechecks its sources.

The formatted copy of a provisional log stays provisional. Its `provisional_manifest.json` keeps the records `provisional_export_eligibility`, `tolerance_resolution` and `warnings` of the canonical run, and its `report.md` repeats the eligibility status and the warnings. The timing recheck of the formatted file uses the tolerance in seconds that the canonical run saved in its `config.json`. No tolerance is converted from a distance at this step.

Before using an exported log for reconstruction:

1. Confirm acquisition identity, chosen signal channels, laser sample selector and the MS first-sample coordinate.
2. For an affine log, use the documented coordinate convention and verify any importer offset or automatic synchronization setting. Do not apply the rate correction twice.
3. For a rate-only log, establish its constant offset and verify that firing-interval endpoints and continuous-trajectory transitions use the same corrected coordinate.
4. Inspect image placement and line transitions against the recorded geometry. Preserve missing data as missing.
5. Save and reopen the receiving workspace. Record the application version, actual settings, imported-file hashes and checks of the resulting time mapping, image and exported values.

The flags `importer_roundtrip_tested`, `gui_workflow_tested` and `receiver_interval_trajectory_agreement_verified` remain false unless the corresponding receiving-application checks were actually performed. Internal export checks must not be described as an application round trip.

The shipped chromium-grating example reproduces timing and corrected-log readback. It contains no saved workspaces. Evidence from a different imported file or prior manuscript reconstruction does not validate newly exported bytes automatically.
