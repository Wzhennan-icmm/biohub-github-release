# Legacy provenance

Original Python, R, AWK, and shell scripts remain source material for migration
and parity validation. Import only source files: exclude virtual environments,
IDE metadata, private inputs, generated figures, and binaries.

For each imported script, record original path, SHA256, owner/license, mapped
BioHub command ID, known behavior changes, and golden test fixture. Third-party
code belongs under `third_party/` with its original license, not under `legacy/`.
