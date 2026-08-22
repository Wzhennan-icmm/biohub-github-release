# AI usage disclosure

BioHub may use generative AI for code scaffolding, refactoring, tests, and
documentation. Before release, maintainers must record tool/model, assisted files,
assistance type, and human validation performed.

AI output is never accepted without human review, domain validation, automated
tests, and author responsibility for correctness. `paper.md` must include this
same disclosure in the journal-required AI usage section.

## v0.4.0 development record (2026-08-14)

- Tool: OpenAI Codex; exact hosted model deployment was not recorded in repository metadata.
- Assistance: recipe architecture, Rust CLI/provenance code, Python/R workflow
  helpers, tests, migration ledger, containers, CI, and documentation.
- Automated validation: Rust format/Clippy/tests; Python helper regressions;
  Python/R syntax checks; all-recipe Snakemake dry-runs; four synthetic low-risk
  workflow executions; checksum-bundle verification.
- Remaining human work: maintainers retain authorship responsibility. Scientific
  interpretations and experimental recipes still require representative data,
  independent comparison, domain review, and approved golden results before
  publication claims.

## v0.4.0 release-hardening record (2026-08-22)

- Assistance: deterministic output hardening, publication golden tests, release
  evidence register, source-entry consolidation, and gated GitHub workflows.
- Human-reviewed source references: sanitized historical scripts already present
  in local migration workspace; no private biological data copied into repository.
- Automated findings corrected: FASTA headers counted as GC sequence, previous
  peptide assigned to next header, invalid GFF Parent reference, non-inclusive
  CDS lengths, dropped first expression row, and infinite second-hit sentinel.
- Release policy: domain, recipe, and visual entries remain explicitly pending;
  automation cannot self-approve them and formal release stays blocked.
