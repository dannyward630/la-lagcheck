# Study summary tables

- `Supplementary_Data_S1.csv`: lag estimates, conditional uncertainty, held-out
  prediction, acquisition duration, scan speed, displacement and diagnostic
  settings/status for all 29 acquisitions (R01 to R29).
- `Supplementary_Data_S2_software_configurations.csv`: paired software outcomes
  for the 37 evaluated configurations with and without the input-placement trial.
- `Supplementary_Data_S3_boundary_checks.csv`: post hoc checks of predictions
  at record-truncated search boundaries.

The repository file `Supplementary_Data_S1.csv` corresponds to the journal
attachment named `Supplementary_Data_S1_R01-R29.csv`; they contain the same
public table. The shorter repository filename is used by the Figure 3 script.

Scientific values and status fields are copied unchanged from the manuscript's
supporting tables. Private source paths, source-file hashes, legacy acquisition
codes and acquisition dates are omitted from the public summary copies. The
software-configuration and boundary tables use anonymous acquisition identifiers.
Among the isotope channels of the biological samples only C13 and Fe56 are named;
a matching signal that includes other channels is given by its number of channels.
In `Supplementary_Data_S1.csv`, the rate fields come from analyses that used the channel sum as the matching signal and every scan line. The `preset_automatic_*` and `preset_declared_*` status fields, the matching signal selected by default (`preset_matching_signal`), its automatic iterated rate and the number of scan lines left out at the start of the record are those of release 0.7.0 under its default settings, with the automatic or declared starting alignment. The laser-log screen flags one line-change interval in each of R05 and R07 (`flagged_line_change_intervals`). The counts of supported acquisitions were obtained after the defaults had been chosen from these records and are retrospective. The S2 and S3 tables record runs of earlier versions (S2 states 0.4.0). Those runs predate the screen: rerun with release 0.6.0, the R05 configuration listed with a verified log in S2 passes the six checks, has one flagged interval and gives a provisional log.
The R01 figure uses the iterated rate of the default 0.7.0 analysis of R01 (Fe56 as the matching signal); the R01 row of the table of the 29 acquisitions keeps the channel-sum, every-line single-pass rate of that table. The two rates differ.

Units are specified in the column names. `ppm` is the software's rate-field name
and is numerically equivalent to microseconds per second. A diagnostic label
describes the stated timing checks, not independent spatial accuracy. These
tables do not contain raw biological signals from acquisitions other than R01.
