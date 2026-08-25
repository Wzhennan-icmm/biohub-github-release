# Statistical recipe validation pack

Public CC0-1.0 synthetic data cover inventory IDs `049`, `050`, `053`, and
`072`. Data contain no human subjects, locations, accessions, or private project
identifiers.

Reference environment is BioHub omics image: R 4.5.3, vegan 2.7-5, and DESeq2
1.50.2. Acceptance tolerates mathematically equivalent floating-point and RDA
axis-sign changes; identifiers, counts, set membership, and permutation results
remain strict.

Build and verify:

```bash
python3 tools/validation_review.py build --pack statistics
python3 tools/validation_review.py verify --pack statistics
```

Human approval must separately confirm model design, contrast, universe, and
multiple-testing interpretation.
