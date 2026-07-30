# BioHub

Rust implementation lives in `biohub-rs`.

## Build and run

```bash
cd biohub-rs
cargo build --release
./target/release/biohub --help
```

Compatibility script examples:

```bash
git clone https://github.com/Wzhennan-icmm/biohub-github-release.git
cd biohub-github-release
python3 run_change_HJJN_geneName.py -i map.txt -o out.txt
python3 run_filter_NCBI_gff_to_get_gene.py --input in.gff3 --output out.gff3
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

## Compatibility entry files

`run_*.py` wrappers and `gff_longest.py` now call the Rust launcher `biohub-rs/run-biohub.sh` and keep old command names usable.

## Legacy content retained

- `dotPlotly`, `*.R`, and `*.awk` are unchanged.
