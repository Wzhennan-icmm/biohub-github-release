# BioHub

BioHub unifies reproducible command-line utilities used in plant-genome assembly,
annotation, comparative genomics, population analysis, and expression workflows.
Core commands are implemented in Rust. R-backed visual and statistical commands
are added through an explicit, versioned backend contract.

中文完整文档：[BioHub v0.4 功能说明书](docs/USER_GUIDE.zh-CN.md)

## Install

### Native core

```bash
git clone https://github.com/Wzhennan-icmm/biohub-github-release.git
cd biohub-github-release/biohub-rs
cargo build --release --locked
./target/release/biohub --help
```

Release archives provide Linux/macOS binaries plus R backends, recipes,
compatibility wrappers, licenses, and metadata. Verify release SHA256 files before use.

### Core container

```bash
docker build -t biohub:local .
docker run --rm biohub:local catalog --format json
docker run --rm -v "$PWD:/work" -w /work biohub:local doctor
```

Core container includes `Rscript`, `samtools`, `mafft`, `pal2nal.pl`, recipes, and
fixtures. Domain recipes need dependencies listed by `biohub doctor --recipe ID`.
Local builds use mutable upstream package repositories and are not bitwise locked;
for publication, retain released image digest plus exact environment export.

## Stable command surface

```bash
biohub catalog
biohub catalog --format json --kind command
biohub catalog --format json --kind recipe
biohub run annotation-vcf --help
biohub run annotation-vcf --reference ref.fa --gff genes.gff3 --vcf calls.vcf --output calls.tsv
biohub doctor
biohub recipe list
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

## Reproducible recipes

Thirteen experimental Snakemake recipes cover comparative genomics, assembly,
population analysis, de novo mutation rates, RNA-seq, enrichment, and microbiome
RDA. Each pack includes config schema/template, workflow, validation, summaries,
input checksums, logs, README, and archive target.

```bash
biohub recipe init selection-branch-site --workdir analysis-config
# Edit config.yaml; copied config.schema.yaml defines accepted fields.
# REQUIRED/null values are placeholders.
biohub doctor --recipe selection-branch-site --strict
biohub recipe validate selection-branch-site --config analysis-config/config.yaml
biohub recipe run selection-branch-site \
  --config analysis-config/config.yaml --workdir runs/selection-001 --cores 8
biohub recipe report --workdir runs/selection-001
```

Runs reject overwrite. `--resume` requires unchanged recipe ID and config SHA256.
`versions.tsv` records every declared dependency as versioned, unavailable, or
available without a safe version probe. `provenance.json` records workflow SHA256
and container hint; set `BIOHUB_CONTAINER_DIGEST` to preserve runtime image digest.
`recipe.sources.sha256` covers packaged workflow, schema, and helper sources.
`checksums.sha256` covers immutable bundle files and excludes mutable `run.json`.
Packaged Slurm profile uses Snakemake executor plugin and conservative one-job
default; copy it and add site-approved resources before cluster use.

New `biohub run <script-id>` protects primary outputs from accidental overwrite.
Pass `--force` to replace an existing primary output. Legacy `scripts run` and
historic wrapper commands keep prior overwrite behavior through v1.x.

`orthofinder-to-pal2nal` validates unique sequence IDs, matched CDS records, CDS
length divisibility by three, and per-group protein/CDS counts before alignment.
It writes `validation_summary.tsv`, runs PAL2NAL with `-nogap`, and returns non-zero
when any orthogroup is skipped or fails.
The historical misspelling `orthofiner-to-pal2nal` remains an alias through v1.x.

## Dependencies and reproducibility

Run `biohub doctor` before workflows using external tools. Commands declare their
backend and dependencies in catalog JSON. R-backed implementations must use
explicit arguments and must not install packages at runtime; see [R backend
contract](r/README.md).

Use deidentified fixtures and approved golden outputs for validation. Never commit
private genomes, sample identifiers, credentials, or unpublished results.

## R commands

`biohub r list` shows R-backed commands. `dotplot` accepts PAF or MUMmer
coordinates. `psmc-plot` accepts merged `Sample/Time/Ne` trajectories and optional
explicit stage intervals. Both write PDF/PNG with base R and reject existing output
unless `--force` is supplied.

## Development

```bash
cd biohub-rs
cargo fmt --check
cargo clippy --locked -- -D warnings
cargo test --locked
cd ..
python3 tools/validate_recipes.py
python3 -m unittest discover -s tools/tests -v
python3 tools/snakemake_smoke.py
```

See the [Chinese functional manual](docs/USER_GUIDE.zh-CN.md),
[CONTRIBUTING.md](CONTRIBUTING.md), [CHANGELOG.md](CHANGELOG.md),
[CITATION.cff](CITATION.cff), [migration ledger](docs/SCRIPT_MIGRATION.md), and
[AI usage disclosure](docs/AI_USAGE.md).

## License and third-party code

BioHub is MIT licensed. Vendored dotPlotly code retains upstream MIT copyright and
notice; see [NOTICE](NOTICE).
