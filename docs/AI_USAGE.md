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
