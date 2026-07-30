# BioHub

Rust implementation lives in `biohub-rs`.

## Build and run

```bash
cd biohub-rs
cargo build --release
./target/release/biohub --help
```

CLI examples:

```bash
git clone https://github.com/Wzhennan-icmm/biohub-github-release.git
cd biohub-github-release
cd biohub-rs
cargo build --release
./target/release/biohub --help
```

Notes:

- `gff/filter-ncbi` 默认输出 `TA-filtered.gff3`（未传 `-o` 时）。
- `gff/filter-gemoma` 默认输出 `gemoma-longest.gff3`（未传 `-o` 时）。

## Command map

- `rename`
  - `hjjn-genes`
  - `scaffolds`
  - `fasta-scaffolds`
- `blast`
  - `reciprocal`
- `gff`
  - `filter-ncbi`
  - `filter-gemoma`
  - `convert-ty1-hjjn`
- `fasta`
  - `longest-transcript`
- `stats`
  - `coverage-ratio`
  - `hic-matrix-reindex`
  - `wgcna-weight`
- `psmc`
  - `merge`

## Legacy content retained

- `dotPlotly`, `*.R`, and `*.awk` are unchanged.
