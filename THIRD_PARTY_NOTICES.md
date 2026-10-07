# Third-party material and release scope

The authored la-lagcheck source is supplied under the existing MIT `LICENSE`.
That license does not relicense outside data or accompanying third-party files.

`tests/fixtures/thermo_qtegra_ldr_5points.csv` is the unmodified five-observation
format fixture from the pewlib project and its contributors, commit
`c77805375a91f96cfe221652cbe94d536eab1011`. It retains its upstream GNU GPL version 3
terms and full license in `tests/fixtures/thermo_pewlib_LICENSE.txt`. The source
URL, original file hash and retrieval provenance are in the adjacent
`thermo_qtegra_ldr_SOURCE.md`. Do not describe the whole release bundle as uniformly
MIT. This fixture is data used by tests; la-lagcheck does not import pewlib or copy its
parser source. Preserve these notices, the fixture's source and full license if
redistributing the fixture.

The pinned upstream [license](https://github.com/djdt/pewlib/blob/c77805375a91f96cfe221652cbe94d536eab1011/LICENSE)
was checked again on 2026-09-28. Its GPLv3 sections 4 and 5 address preserving
notices with verbatim source copies and distinguish an aggregate of separate
works. This records the files' actual stated terms; it is not a new determination
of institutional ownership or a claim of blanket clearance to publish the paper's
experimental data.

The synthetic generators in `examples/` do not read biological acquisitions.
The optional R005 raw fixture and complete iolite workspaces are absent. A separate,
authorised R01 example contains only anonymous C13 and Fe56 signal values and firing
records, as described in `examples/R01/README.md`. No test or example uploads data,
sends an external message or publishes files.

## Chromium-grating raw data

The two CSV files in `examples/chromium_grating/inputs/` provide the grating acquisition used for reproduction. The MS instrument-path header has been replaced by a neutral description; numerical records and the laser log are unchanged. The software MIT license does not license these experimental data. No additional data license is asserted by this package. See `DATA_SCOPE.md` for the biological-example scope and exclusions.
