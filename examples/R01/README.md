# R01: a representative biological imaging acquisition

This example reproduces the position and local-orientation comparisons in
Figure 3 of *Measuring and correcting time-dependent lag in long LA-ICP-MS
imaging acquisitions*. R01 lasted 9.184 h. Only the authorised C13 and Fe56
signals are distributed, under the anonymous identifier R01. Fe56 was used for
the figure and the post hoc rate analysis.

## Files and provenance

- `R01_signals.csv` contains 505016 stored signal observations, with columns
  `time_s`, `C13` and `Fe56`. Signals are in counts per second. These values were
  extracted from the saved iolite workspace used for Figure 3, without changing
  intensity values. They are not represented as an original vendor-export file.
  The last observation is missing in both channels and remains `nan` in this file.
- `R01_workspace_intervals.csv` contains all 904 saved firing intervals and their
  endpoint coordinates in micrometres. Its times share the signal file's
  anonymous relative coordinate. The workspace's input placement is preserved.
- `R01_native_laser_log.csv` contains only the selected acquisition's original
  5426 laser-log records. Sample labels have been replaced by R01, and timestamps
  have been converted to seconds relative to the first selected record. Other
  cells and record order are unchanged. This separate native representation has
  903 positive-duration intervals and one exactly zero-duration On/Off pair.
  The saved workspace represents the corresponding very short interval as
  approximately 0.48 microseconds. The two representations are kept separate.
- `reference_mapping.json` and `reference_blocks.csv` specify the post hoc mapping
  used for Figure 3. `expected_figure3.json` records the reported comparison.
- `reconstruction.py` and `feature_orientation.py` reproduce the image arrays and
  local-feature extraction. `reproduce_figure3.py` writes the complete comparison.
- `reestimate_posthoc.py` repeats the preset Fe56 diagnosis and the separate
  post hoc boundary treatment, retaining their different outcomes.

Absolute dates, original sample labels, workspace identifiers, image metadata and
local file paths are omitted. A shared integer time-origin shift leaves relative
times unchanged. For exact numerical reproduction, the scripts add a synthetic
power-of-two origin that preserves the original floating-point resolution. This
synthetic value is not an acquisition timestamp. The original saved workspace is
not distributed.

## Reproduce the comparison

Install the package's reference dependencies, then run from its root:

```text
python -B -X utf8 examples/R01/reproduce_figure3.py --output ../R01_figure3
python -B -X utf8 examples/R01/reestimate_posthoc.py --output ../R01_posthoc
```

Use a new output directory for each command. The first command writes PNG, PDF
and SVG figures and a numerical audit. It reconstructs both images from the
signals and firing intervals; no precomputed biological image or mask is included.

The second command excludes the one terminal missing observation from rate
estimation, leaving a single contiguous segment of 505015 observations. No
missing value is replaced or interpolated for estimation. All 904 saved intervals
are retained. The Fe56 rate with the documented post hoc boundary rule is
168.99125403067967 microseconds per second (approximately 169.0 in the paper).
The preset Fe56 diagnosis fails checks, whereas the post hoc analysis is supported
under its stated conditions. The historical channel selection among all saved
channels is not repeated because only C13 and Fe56 are distributed.

The common 15 micrometre reconstruction grid reproduces the two 902 by 1114
arrays used for Figure 3. Expected early and late along-scan shifts are
-0.7907883197673672 and +0.7740541140814957 mm. The selected local feature has
principal-axis angles of 24.993774508375427 and 12.652875688748152 degrees in the
constant-offset and mapped reconstructions, respectively. The difference is
local to the selected feature; it is not a global shear angle or an independently
known anatomical direction.

The all-acquisition panel uses `study_summaries/Supplementary_Data_S1.csv`.
Its R01 single-pass rate differs from the post hoc Fe56 mapping and must not be
substituted for the latter. The summaries retain all 29 acquisitions; the raw
signals of the other biological acquisitions are not included.

## Interpretation and terms

This example reproduces an analysis using post hoc settings. It does not provide
independent spatial validation or change the standard LA-LagCheck diagnostic.
The native event log is supplied for provenance and input-format inspection;
the Figure 3 commands use the saved-workspace representation explicitly.

Authored scripts use the package's existing MIT licence. Experimental data are
separate from that software licence; no additional data-reuse licence is asserted
here. No clinical attributes, other biological channels, complete workspaces,
external-laboratory datasets or AI prompt records are included.
