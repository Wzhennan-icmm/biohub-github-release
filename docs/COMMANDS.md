# Command support policy

`biohub catalog --format json` is machine-readable source of command discovery.
Every record reports ID, original source, status, category, backend, and external
dependencies.

Support levels:

- **stable**: documented inputs/outputs, regression fixture, CI execution.
- **optional**: requires named external executable or R package; `doctor` reports it.
- **example**: educational or demonstration workflow; not scientific-analysis API.
- **deprecated**: legacy wrapper retained through v1.x with migration target.

Current Rust catalog records use legacy status while fixtures are being expanded.
Do not relabel an implementation as stable until its parity fixture and docs exist.
