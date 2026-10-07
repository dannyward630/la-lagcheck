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
The R01 figure uses a post hoc Fe56 iterated rate; the R01 row in the cohort summary
uses the separately documented default-signal single-pass rate.

Units are specified in the column names. `ppm` is the software's rate-field name
and is numerically equivalent to microseconds per second. A diagnostic label
describes the stated timing checks, not independent spatial accuracy. These
tables do not contain raw biological signals from acquisitions other than R01.
