# Domain containers

Core `Dockerfile` stays small. `Dockerfile.domain` builds five domain images from
direct dependency constraints:

```bash
docker build -f containers/Dockerfile.domain \
  --build-arg ENVIRONMENT=comparative \
  -t biohub-comparative:0.4.0 .
```

| Environment | Recipes |
| --- | --- |
| comparative | orthology/codon, CAFE, PAML selection/dating, synteny/SV |
| assembly | T2T assembly evidence |
| population | SNP GWAS and population selection |
| variant | family de novo rate |
| omics | RNA-seq, enrichment, microbiome RDA |

Direct versions were checked 2026-08-14 against package-provider records:
[Snakemake](https://anaconda.org/bioconda/snakemake),
[Slurm executor](https://anaconda.org/bioconda/snakemake-executor-plugin-slurm),
[CAFE](https://anaconda.org/bioconda/cafe),
[PAML](https://anaconda.org/bioconda/paml),
[MAFFT](https://anaconda.org/bioconda/mafft),
[PAL2NAL](https://anaconda.org/bioconda/pal2nal),
[minimap2](https://anaconda.org/bioconda/minimap2),
[SyRI](https://anaconda.org/bioconda/syri),
[PLINK2](https://anaconda.org/bioconda/plink2),
[VCFtools](https://anaconda.org/bioconda/vcftools),
[DESeq2](https://anaconda.org/bioconda/bioconductor-deseq2), and
[vegan](https://anaconda.org/conda-forge/r-vegan).

These YAML files pin direct tools, not complete transitive builds. Published image
digest plus `micromamba list --explicit` export is release lock. Release checklist
requires preserving both. Never treat mutable image tag alone as provenance.
Domain-container workflow uploads explicit export, local image ID, and—on tags—
published registry digest as per-environment metadata artifacts.

`kmer-gwas` is excluded from published images because its external analysis script
requires end-of-life Python 2. Supply an audited isolated legacy environment and
record its digest; do not add Python 2 to other domain images.

Slurm profile also needs cluster-provided `sbatch` and site-approved resource
settings. Image does not infer account or partition.
