# Changelog

All notable changes follow [Semantic Versioning](https://semver.org/).

## 0.3.0 - Unreleased

- Add stable `catalog`, `run`, and `doctor` command surfaces.
- Add JSON command catalog metadata and external-dependency preflight.
- Add release, container, contribution, and citation scaffolding.
- Validate protein/CDS parity and codon frames before MAFFT/PAL2NAL conversion.
- Write per-orthogroup validation status and use gap-free PAL2NAL output.
- Package R backends and compatibility wrappers in native release archives.
- Exercise base-R dotplot generation in container CI.
- Route R-backed dotplot through unified `biohub run` interface.
- Preserve nested command help instead of replacing it with global help.
- Parse whitespace- and pipe-delimited MUMmer coordinate tables.
- Reject release tags that do not match Cargo package version.
- Add process-level CLI tests for help exit codes, JSON catalog metadata, and unknown commands.
- Add a comprehensive Chinese functional manual covering all 54 catalog commands.
- Add and exercise a synthetic MUMmer coords dotplot fixture.

## 0.2.1

- Rust implementation of legacy BioHub command catalog.
