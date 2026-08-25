# Orthology and codon validation pack

Public CC0-1.0 synthetic three-taxon data cover inventory IDs `030`, `033`,
`034`, and recipe `043`. Variable synonymous codons (`GCT/GCC/GCA`) distinguish
the two fourfold-degeneracy contracts. Invariant `GGT` supplies strict sites.

Pinned release reference environment:

- MAFFT 7.525: `mafft --maxiterate 1000 --localpair`
- PAL2NAL 14.1: `pal2nal.pl protein.aln cds.fa -output paml -nogap`
- BioHub recipe config: `configs/comparative-orthology-codon.json`

Build and verify:

```bash
python3 tools/validation_review.py build --pack orthology-codon
python3 tools/validation_review.py verify --pack orthology-codon
```

Verification translates each codon alignment independently with standard genetic
code. Human reviewer must still approve reciprocal-set semantics and alignment.
