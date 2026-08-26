# Command support policy

`biohub catalog --format json` is machine-readable source of command and recipe
discovery. Every record reports ID, kind, status, domain, backend, dependencies,
version, and license. Command records retain original source; recipe records add
workflow, schema, and container.

User-facing syntax, inputs, outputs, limitations, and tutorials are documented in
the [Chinese functional manual](USER_GUIDE.zh-CN.md).

Support levels:

- **stable**: documented inputs/outputs, regression fixture, CI execution.
- **optional**: requires named external executable or R package; `doctor` reports it.
- **example**: educational or demonstration workflow; not scientific-analysis API.
- **deprecated**: legacy wrapper retained through v1.x with migration target.
- **experimental recipe**: reproducible workflow contract exists, but domain-approved golden outputs remain pending.

Current Rust catalog records use legacy status while fixtures are being expanded.
Do not relabel an implementation as stable until its parity fixture and docs exist.
Migration decisions and deferred scripts are tracked in
[`SCRIPT_MIGRATION_MATRIX.tsv`](SCRIPT_MIGRATION_MATRIX.tsv).
