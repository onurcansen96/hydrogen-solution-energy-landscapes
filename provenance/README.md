# Provenance and Publication Checks

The scientific records were supplied by the dataset author. Publication
packaging must not be confused with running new DFT calculations.

- `author_source_manifest.json` is the original source manifest, unmodified.
  Its README/converter hashes describe the author's export before the public
  documentation edits. Its statement about absent volumetric files predates
  inclusion of the CHGCAR release assets. Use the current root inventory for
  public release contents.
- `author_validation_report.json` is the original report, unmodified. It
  includes checks against external files that are not distributed here.
- `release_audit.json` records publication-time byte comparisons against the
  author's inputs and scientific JSON, plus the isolated converter smoke test.
  These are integrity/consistency checks, not independent validation of DFT
  accuracy or of all manuscript figures.

The raw 25 dumps and descriptor table, all 50 per-structure JSON files,
combined conceptual dataset, and computational settings are preserved
byte-for-byte. Scientific property values are not changed.

Publication edits are limited to documentation/licenses/examples, the
charge-density inventory, clarified CD resource wording, and updated
packaging manifests/checksums. The converter's scientific calculations are
unchanged; its maintenance update describes release assets and requires
explicit `--overwrite` before replacing existing generated output.

`generated_metadata/validation_report.json` retains the original author
report and timestamp. `generated_metadata/source_manifest.json` records a
separate `publication_packaging` entry and refreshes hashes for edited support
files. `generated_metadata/manifest_index.json` refreshes output hashes while
retaining the original scientific-generation timestamp and version.

The unpublished manuscript, representative INCAR, and selected-descriptor
cross-check CSV are historical external references, not downloadable contents
of this release. Their original hashes identify them without exposing local
computer paths. No citation, DOI, or missing physical result has been invented.
