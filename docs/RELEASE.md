# Release checklist

1. Run `python3 tools/validate_release.py`. Resolve every inconsistency in
   migration matrix and validation register.
2. Obtain named human approval for every domain, recipe, and visual `pending`
   row. Never convert synthetic execution into human approval.
3. Update matrix and `validation/reviews.tsv` in same reviewed commit. Formal
   command `python3 tools/validate_release.py --release --tag vX.Y.Z` must pass.
4. Run Rust format, Clippy, all tests, external samtools goldens, recipe static
   validation, Python regressions, all 13 DAG dry-runs, and four safe executions.
5. Merge only green reviewed commit to `main`; manually dispatch Domain
   containers with `publish=false` and release version as image tag.
6. Update version, release date in `CHANGELOG.md`, `CITATION.cff`, documentation,
   compatibility notes, tested platforms, and known limitations.
7. Create annotated immutable `vX.Y.Z` tag. Release workflow re-runs formal gate,
   builds three native archives, validates their SHA256/content, builds five
   domain images, and publishes GitHub Release only after all jobs pass.
8. Download every release asset, verify SHA256, run `--version`, `catalog`,
   `doctor`, example plotting, and one safe Recipe on a clean machine.
9. Record five immutable GHCR digests and exact environment exports; never
   advertise mutable-only image tags.
10. Archive tagged release with Zenodo and add generated DOI to `CITATION.cff` in
    next metadata release. For JOSS, verify public-history timing, paper metadata,
    impact citations, and AI disclosure.
