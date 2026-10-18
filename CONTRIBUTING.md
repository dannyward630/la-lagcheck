# Contributing to LA-LagCheck

LA-LagCheck is maintained by Linhan (Leo) Li at the MATFab Core Facility, University of Iowa. Bug reports, questions and suggestions are welcome, including from people who only want to try the program on their own data.

## Reporting a problem

Open an issue at <https://github.com/yzoe236/la-lagcheck/issues> and include:

- the release you downloaded or the commit you used (the `version` field of the report also shows it);
- the command you ran and its full console output;
- the laser and ICP-MS models and the versions of their acquisition software;
- whether you can share a small example file.

Please do not attach data you are not allowed to share. That includes anything that identifies a patient, an animal or a collaborator's sample. A short description of the file layout is often enough.

## Proposing a change

1. For anything larger than a typo, open an issue first so we can agree on the approach.
2. Fork the repository and work on a branch.
3. Install the tested environment from the repository root:
   `python -m pip install -r requirements-reference.txt`
4. Run the tests from the repository root:
   `python -B -X utf8 -m pytest -q -p no:cacheprovider tests`
5. Open a pull request that explains what changed and why, and add a line to `CHANGELOG.md`.

A change to the timing model or to a default behaviour needs a test that shows the old and the new result on the same input.

## Support and maintenance

- The maintainer aims to reply to new issues within two weeks.
- Releases are tagged on GitHub with a changelog entry. Selected releases are archived on Zenodo; the concept DOI <https://doi.org/10.5281/zenodo.23324268> always resolves to the latest archived release.
- `FILES.json` and `scripts/verify_package.py` describe the files of the latest tagged release. Files on the `main` branch can change between releases; use a tagged release when you need a fixed set of files.
- Automated tests run on every push to `main` and on every pull request.
