# BioHub

BioHub unifies reproducible command-line utilities used in plant-genome assembly,
annotation, comparative genomics, population analysis, and expression workflows.
Core commands are implemented in Rust. R-backed visual and statistical commands
are added through an explicit, versioned backend contract.

中文完整文档：[BioHub v0.3 功能说明书](docs/USER_GUIDE.zh-CN.md)

## Install

### Native core

```bash
git clone https://github.com/Wzhennan-icmm/biohub-github-release.git
cd biohub-github-release/biohub-rs
cargo build --release --locked
./target/release/biohub --help
```

Release archives provide Linux/macOS binaries plus R backends, compatibility
wrappers, licenses, and metadata. Verify release SHA256 files before use.

### Reproducible container

```bash
docker build -t biohub:local .
docker run --rm biohub:local catalog --format json
docker run --rm -v "$PWD:/work" -w /work biohub:local doctor
```

Container includes `Rscript`, `samtools`, `mafft`, and `pal2nal.pl`. Hamstr is an
optional host-provided dependency.

## Stable command surface

```bash
biohub catalog
biohub catalog --format json
biohub run annotation-vcf --help
biohub run annotation-vcf --reference ref.fa --gff genes.gff3 --vcf calls.vcf --output calls.tsv
biohub doctor
biohub run dotplot --input alignments.paf --output dotplot.pdf --format paf
```

`biohub scripts catalog` and `biohub scripts run <script-id>` remain supported
through v1.x. Existing direct command groups remain available:

- `rename`: HJJN genes, scaffold IDs, FASTA scaffold IDs
- `blast`: reciprocal BLAST filters
- `gff`: NCBI/GeMoMa filtering and coordinate conversion
- `fasta`: longest transcript selection
- `stats`: coverage ratio, Hi-C matrix reindexing, WGCNA weights
- `psmc`: merge plotting inputs

Run `biohub catalog` for complete IDs, source provenance, status, and dependencies.

New `biohub run <script-id>` protects primary outputs from accidental overwrite.
Pass `--force` to replace an existing primary output. Legacy `scripts run` and
historic wrapper commands keep prior overwrite behavior through v1.x.

`orthofiner-to-pal2nal` validates unique sequence IDs, matched CDS records, CDS
length divisibility by three, and per-group protein/CDS counts before alignment.
It writes `validation_summary.tsv`, runs PAL2NAL with `-nogap`, and returns non-zero
when any orthogroup is skipped or fails.

## Dependencies and reproducibility

Run `biohub doctor` before workflows using external tools. Commands declare their
backend and dependencies in catalog JSON. R-backed implementations must use
explicit arguments and must not install packages at runtime; see [R backend
contract](r/README.md).

Use deidentified fixtures and approved golden outputs for validation. Never commit
private genomes, sample identifiers, credentials, or unpublished results.

## R commands

`biohub r list` shows R-backed commands. `biohub run dotplot` and
`biohub r run dotplot` are equivalent. `dotplot` accepts PAF or MUMmer
coordinate rows and writes PDF/PNG with base R only. It rejects existing outputs
unless `--force` is supplied.

## Development

```bash
cd biohub-rs
cargo fmt --check
cargo clippy --locked -- -D warnings
cargo test --locked
```

See the [Chinese functional manual](docs/USER_GUIDE.zh-CN.md),
[CONTRIBUTING.md](CONTRIBUTING.md), [CHANGELOG.md](CHANGELOG.md),
[CITATION.cff](CITATION.cff), and [AI usage disclosure](docs/AI_USAGE.md).

## License and third-party code

BioHub is MIT licensed. Vendored dotPlotly code retains upstream MIT copyright and
notice; see [NOTICE](NOTICE).
