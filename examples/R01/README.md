# R01: a representative biological imaging acquisition

This example reproduces the position and local-orientation comparisons in
Figure 3 of *Measuring and correcting time-dependent lag in long LA-ICP-MS
imaging acquisitions*. R01 lasted 9.184 h. Only the authorised C13 and Fe56
signals are distributed, under the anonymous identifier R01. The figure uses the
mapping from the default LA-LagCheck 0.7.0 diagnosis, in which Fe56 is the
matching signal.

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
- `reference_mapping.json` and `reference_blocks.csv` specify the mapping used for
  Figure 3: the diagnosis with the default settings and automatic starting
  alignment. Its tolerance, which has no default and must be given, was set to
  0.2 s. `expected_figure3.json` records the reported comparison.
- `reconstruction.py` and `feature_orientation.py` reproduce the image arrays and
  local-feature extraction. `reproduce_figure3.py` writes the complete comparison
  and checks its numbers against `expected_figure3.json`.
- `reestimate_default.py` repeats the default diagnosis from the distributed
  signals and checks it against `reference_mapping.json` and
  `reference_blocks.csv`.

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
python -B -X utf8 examples/R01/reestimate_default.py --output ../R01_default
```

Use a new output directory for each command. The first command writes PNG, PDF
and SVG figures and a numerical audit. It reconstructs both images from the
signals and firing intervals; no precomputed biological image or mask is included.

The figure is drawn with Arial. If Arial is not installed, the script uses DejaVu
Sans, which is distributed with matplotlib, and prints one line that says so.
The numerical audit is the same in both cases. Without Arial the lettering
differs from the published figure, so the PNG is not pixel-identical to it.

The second command excludes the one terminal missing observation from rate
estimation, leaving a single contiguous segment of 505015 observations. No
missing value is replaced or interpolated for estimation. The matching signal is
chosen by the default rule: the C13+Fe56 sum, C13 and Fe56 are compared by the
peak of the automatic starting-alignment search (0.428, 0.441 and 0.484), and
Fe56 is used. With every scan line, the first-block correlation peak lies at the
search boundary set by the start of the signal record in the 12- and 15-block
alternatives. The default rule `skip_record_edge_lines` therefore repeats the
diagnosis without the three scan lines that start within 30 s of the first
recorded sample. All fourteen checks then pass, the status is `supported`, all
ten blocks are retained, and the iterated rate is 169.05931388632462
microseconds per second (approximately 169.1 in the paper). The single-pass rate
of the same run is 169.7961267682242. With the saved input placement declared
(`alignment_offset_s = 0`), the status is also `supported` and the iterated rate
is 169.05999859650836. The mapping is applied to all 904 saved intervals. The
selection among all saved channels is not repeated because only C13 and Fe56 are
distributed.

The option `exclude_record_boundary_blocks` still exists and remains off by
default. With `skip_record_edge_lines` switched off, it leaves the first block
out of the 12- and 15-block alternatives and gives the status
`supported_excluding_record_boundary_blocks` with a rate of 168.99125403067967;
that analysis was used for an earlier version of the figure.

The common 15 micrometre reconstruction grid reproduces the two 902 by 1114
arrays used for Figure 3. Expected early and late along-scan shifts are
-0.7911068033433237 and +0.7743658580919655 mm. The selected local feature has
principal-axis angles of 24.850215543658283 and 12.590642375382895 degrees in the
constant-offset and mapped reconstructions, respectively. The difference is
local to the selected feature; it is not a global shear angle or an independently
known direction in the sample.

The all-acquisition panel uses `study_summaries/Supplementary_Data_S1.csv`.
Its R01 single-pass rate (171.758550802475, from the archived analysis) differs
from the mapping used here and must not be substituted for it. The summaries retain all 29 acquisitions; the raw
signals of the other biological acquisitions are not included.

## Interpretation and terms

This example reproduces an analysis with the default settings. The two defaults
involved were chosen after the outcomes of the study acquisitions had been
examined. The example does not provide independent spatial validation.
The native event log is supplied for provenance and input-format inspection;
the Figure 3 commands use the saved-workspace representation explicitly.

Authored scripts use the package's existing MIT licence. Experimental data are
separate from that software licence; no additional data-reuse licence is asserted
here. No sample attributes, other channels, complete workspaces,
other datasets or AI prompt records are included.
