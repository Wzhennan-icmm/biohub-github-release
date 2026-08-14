# Contributing

## Before opening a pull request

1. Open an issue for behavior changes or new commands.
2. Preserve existing command IDs and legacy aliases through v1.x.
3. Add command-registry metadata, documentation, and regression fixtures.
4. Run `cargo fmt --check`, `cargo clippy --locked -- -D warnings`, and `cargo test --locked`.
5. For recipes, run `python3 tools/validate_recipes.py`, helper regressions, and
   `python3 tools/snakemake_smoke.py`; record any required domain execution.

## Command acceptance rules

Every stable command needs explicit inputs, outputs, failure behavior, dependency
metadata, a small deidentified fixture, and a golden output approved by a domain
reviewer. Do not add automatic package installation, hard-coded absolute paths,
or silently overwritten output files.

Experimental recipes additionally need schema/template parity, explicit scientific
thresholds, dependency/container metadata, input checksums, logs, run provenance,
archive target, synthetic DAG fixture, and a migration-ledger update when replacing
an existing script. Promotion from experimental needs representative deidentified
data and domain-approved golden results.

## Reporting problems

Use GitHub Issues with command ID, BioHub version, operating system, exact
command, sanitized input description, expected result, and actual stderr.
