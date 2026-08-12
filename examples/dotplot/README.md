# Dotplot smoke-test fixture

```bash
biohub r run dotplot --input input.paf --output dotplot.pdf --format paf
biohub run dotplot --input input.coords --output dotplot.png --format coords
```

Fixture contains three synthetic PAF alignments covering forward and reverse
strands and three synthetic MUMmer `show-coords` rows. Successful runs write
non-empty PDF or PNG output.
