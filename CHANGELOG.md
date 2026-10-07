# LA-LagCheck 0.5.0

This local research-preview release adds explicit New Wave event-log parsing, conditional analysis for short line or spot records, and coordinate/timing inspection views. See NEW_WAVE.md for input selection, analysis conditions and reproducible commands. The original long-acquisition estimator and its reference thresholds remain available.

- Sequence/sample identity and original event record numbers are preserved. Repeated On events are audited as state updates. Spatial views use stated first-On/Off endpoint assumptions; stale intermediate coordinates and inherited positioning speeds are not treated as measured ablation trajectories or scan speeds.
- Short records report available, failed and unexecuted checks separately. Unverified automatic peak candidates are labelled as such, including when no rate can be fitted. A failed diagnostic does not suppress readable inspection views. Conditional corrected logs remain provisional.
- Iolite-compatible export can omit only zero-duration event pairs whose prior exclusion is tied to unchanged original input and provenance. It records every omitted row. It can materialize inherited sample/sequence labels before sorting while checking event membership, pairing and all other non-time cells. This does not claim actual receiving-application acceptance.

Validation includes the public otolith and gabbro examples, internal fixed-mapping readback of provisional logs when a rate exists, a chromium-grating control, and software tests. A 69-line sample-only otolith selection did not provide four usable blocks and therefore received inspection views without an estimated correction. Physical time-origin and spatial accuracy remain separate from successful parsing or plotting.

Candidate validation reports may display version 0.4.0 because the version metadata was changed only after the functional code was frozen. Project release records bind the unchanged functional source hashes to this 0.5.0 copy and record a fresh full-suite validation. No public repository upload or external publication is implied by this local directory.
