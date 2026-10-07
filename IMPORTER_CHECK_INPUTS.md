# Receiving-application checks

`iolite-export` formats a verified canonical laser log as an ESL CSV and checks millisecond timestamps, stable event ordering, firing-event pairing, unchanged non-time fields and the declared fixed time mapping. It does not open iolite or HDIP, save a reconstruction workspace, or validate an image in either application.

Keep the original inputs, configuration, canonical export, receipts and SHA256 hashes together. The original input paths must remain accessible when the ESL formatting step rechecks its sources.

Before using an exported log for reconstruction:

1. Confirm acquisition identity, chosen signal channels, laser sample selector and the MS first-sample coordinate.
2. For an affine log, use the documented coordinate convention and verify any importer offset or automatic synchronization setting. Do not apply the rate correction twice.
3. For a rate-only log, establish its constant offset and verify that firing-interval endpoints and continuous-trajectory transitions use the same corrected coordinate.
4. Inspect image placement and line transitions against the recorded geometry. Preserve missing data as missing.
5. Save and reopen the receiving workspace. Record the application version, actual settings, imported-file hashes and checks of the resulting time mapping, image and exported values.

The flags `importer_roundtrip_tested`, `gui_workflow_tested` and `receiver_interval_trajectory_agreement_verified` remain false unless the corresponding receiving-application checks were actually performed. Internal export checks must not be described as an application round trip.

The shipped chromium-grating example reproduces timing and corrected-log readback. It contains no collaborator workspaces. Evidence from a different imported file or prior manuscript reconstruction does not validate newly exported bytes automatically.
