# Annotation and coordinate validation pack

Public synthetic GFF3 records exercise BioHub inventory IDs `014`, `015`, and
`016`. Fixture license is CC0-1.0. No private path, species sample, or project
identifier is present.

Reference convention:

- GFF3 coordinates are 1-based closed intervals.
- Forward transform: `new = old - scaffold_start + chromosome_start`.
- Reverse transform: `new_start = scaffold_end - old_end + chromosome_start`;
  `new_end = scaffold_end - old_start + chromosome_start`; strand flips.
- Unmapped or out-of-range rows must be preserved or recorded in explicit audit
  output. Silent deletion fails review.

Build and independently recheck byte-exact reference outputs:

```bash
python3 tools/validation_review.py build --pack annotation-coordinates
python3 tools/validation_review.py verify --pack annotation-coordinates
```

Human approval remains separate. Inspect `review.md`, then explicitly report
`批准 annotation-coordinates` or reject with differences.
