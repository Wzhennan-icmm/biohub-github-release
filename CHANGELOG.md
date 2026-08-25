# Changelog

All notable changes follow [Semantic Versioning](https://semver.org/).

## 0.4.0 - Unreleased

- Add 13 schema-validated Snakemake recipe packs for comparative genomics,
  assembly, population genomics, de novo mutation rates, RNA-seq, enrichment,
  and microbiome RDA.
- Add `biohub recipe list|describe|init|validate|run|report`, local/Slurm profiles,
  config-hash resume protection, run state, provenance, checksums, and archives.
- Unify command and recipe discovery in catalog JSON with kind, domain, version,
  dependency, container, and license metadata.
- Add recipe-specific dependency preflight, including R package checks.
- Add canonical `orthofinder-to-pal2nal` ID while retaining misspelled alias.
- Integrate FASTA-length-aware GFF3 filtering with model propagation.
- Add parameterized base-R PSMC trajectory plotting without hard-coded stages.
- Add sanitized 91-item script migration ledger and private inventory tool.
- Package recipes in native archives and core container; add static recipe CI.
- Correct assembly coverage denominators to include unaligned FASTA records; add
  strict PAF length/coordinate checks.
- Add regression coverage for SyRI annotation columns, whitespace-delimited
  MCMCTree chains, PLINK2 hybrid results, and windowed nucleotide diversity.
- Expand run provenance with declared dependency versions, workflow checksum,
  container digest hook, stable bundle checksums, and managed-path collision checks.
- Build five pinned domain images, preflight every included recipe dependency,
  and execute four synthetic low-risk workflows in image CI.
- Add a publication-oriented Chinese recipe manual covering every config field,
  input contract, output, recovery path, and reporting checklist.
- Add deterministic publication goldens for migrated commands, external BAM
  fixtures, visual fingerprints, and a non-bypassable release evidence gate.
- Fix FASTA GC header counting, longest-peptide record attribution, GFF Parent
  references, inclusive CDS lengths, first-row expression loss, and best-hit
  sentinel handling.
- Make map-backed outputs deterministic and reject stale local launcher binaries.
- Gate GitHub Release on validation readiness, three native archives, and five
  tested immutable domain containers.
- Add four public CC0 validation packs with explicit numeric, codon, coordinate,
  and visual acceptance contracts; evidence generation cannot self-approve.
- Prevent SVG scatter tick and chromosome-label clipping with adaptive numeric
  formatting and dedicated custom-axis margins.
- Add structured CFF author/affiliation/ORCID metadata and isolated schema
  validation without downgrading Snakemake's JSON Schema runtime.
- Split release automation into non-publishing manual preflight and annotated-tag
  publication paths; enforce tag target equals current `main`.

## 0.3.0 - Development snapshot

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
- Add a comprehensive Chinese functional manual covering the command catalog.
- Add and exercise a synthetic MUMmer coords dotplot fixture.

## 0.2.1

- Rust implementation of legacy BioHub command catalog.
