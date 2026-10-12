# Experimental data and release scope

This release contains the chromium-grating example and one representative
biological acquisition, R01. The R01 numerical data are restricted to C13 and
Fe56 and anonymised firing records. Acquisition summaries remain available for
all 29 routine acquisitions, together with software-configuration and boundary
check tables. Other biological raw data, complete iolite workspaces, identities,
sample metadata and AI prompt records are excluded.

The grating signal file retains its numerical observations and acquisition
metadata except for the local instrument-path header, which has been replaced
by a neutral description. The R01 signal file is a numerical extraction from
the saved workspace used for the paper, not an original vendor export. Its
README explains the time-origin removal, missing final value and distinction
between native and saved-workspace laser intervals.

The existing MIT licence covers authored software and documentation. It does
not automatically grant a licence to the experimental data or third-party test
material. No additional experimental-data reuse licence is asserted by this
release. The separately attributed Thermo format fixture retains its GPLv3
licence and source notice in `THIRD_PARTY_NOTICES.md`.

The manuscript is not yet represented here as a published article or assigned a
DOI. Repository availability and numerical reproduction are separate from journal
acceptance and from independent verification of image geometry.
