# PSMC plot fixture

Synthetic, non-biological values for CLI and container smoke tests:

```bash
biohub run psmc-plot --input merged.tsv --output plot.pdf \
  --x-scale log10 --y-scale log10 --stages stages.tsv
```
