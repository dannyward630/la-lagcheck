# Thermo Qtegra LDR fixture provenance

`thermo_qtegra_ldr_5points.csv` is an unmodified copy of the five-observation
test fixture in the original pewlib developer repository, `djdt/pewlib`.
Credit: the pewlib project and its contributors. This is a developer-supplied
format example, not a Thermo Fisher supplied validation dataset.

- Repository: https://github.com/djdt/pewlib
- Commit: `c77805375a91f96cfe221652cbe94d536eab1011`
- Source: https://github.com/djdt/pewlib/blob/c77805375a91f96cfe221652cbe94d536eab1011/tests/data/csv/thermo_icap_ldr/test_icap_ldr_1.csv
- Raw bytes: https://raw.githubusercontent.com/djdt/pewlib/c77805375a91f96cfe221652cbe94d536eab1011/tests/data/csv/thermo_icap_ldr/test_icap_ldr_1.csv
- Retrieved: 2026-09-19
- SHA256: `4878036640e2e15d537e2b9e73b550e161ec40a6d31a3b8df2cfd50a3f7da722`
- Size: 3270 bytes; original CRLF line endings retained.
- Upstream license: GNU General Public License version 3; the unmodified
  full upstream LICENSE is retained as `thermo_pewlib_LICENSE.txt`.
- License source: https://github.com/djdt/pewlib/blob/c77805375a91f96cfe221652cbe94d536eab1011/LICENSE

The file identifies Qtegra 2.7.2425.65, File Version 1 and iCAP RQ. It has a
single `Time` column, `31P`, `153Eu`, and `182W` signal columns, a dwell/xcal
metadata row, and a trailing empty field. Tests explicitly configure
`time_unit="s"`; the `Time` header does not state its unit. The fixture tests
file recognition and numerical parsing only. Five observations cannot
validate drift estimation, long acquisitions, other software versions,
other instrument models, or downstream iolite reconstruction.

The LDR adapter was implemented from the observed file structure, without
copying pewlib parser source code. Ordinary Qtegra rows/columns exports can
have separate times for each isotope and are outside this adapter's scope.
