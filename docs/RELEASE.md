# Release checklist

1. Build and verify all public representative packs:
   `python3 tools/validation_review.py build --pack all`, then
   `python3 tools/validation_review.py verify --pack all`.
2. Obtain explicit named human approval for each pack. Never convert automated
   execution into human approval. Save signed records, then update migration
   matrix and `validation/reviews.tsv` in same reviewed commit.
3. Replace `Unreleased` with actual release date in `CHANGELOG.md`; add same
   `date-released` to `CITATION.cff`. Do not predict release date.
4. Run `python3 tools/validate_release.py --release --tag vX.Y.Z`, isolated
   `cffconvert --validate`, Rust format, Clippy, all tests, external samtools
   goldens, recipe static validation, all 13 DAG dry-runs, and four safe runs.
5. Merge only green reviewed PR to `main` using merge commit. Record exact `main`
   SHA; do not tag feature branch or unreviewed commit.
6. From `main`, manually dispatch **Release** workflow with `release_tag=vX.Y.Z`.
   Dispatch is preflight-only: it builds and verifies three native archives and
   five domain images, uploads artifacts/environment exports, and publishes
   neither GHCR images nor GitHub Release.
7. Inspect preflight summary and artifacts. Verify three archives/SHA256 files,
   five environment exports, and five local image IDs.
8. Create annotated, non-cryptographically-signed `vX.Y.Z` tag on exact recorded
   `main` SHA and push tag. Workflow rejects lightweight tags or tags not pointing
   to current `main`.
9. Tag workflow repeats formal gate and builds. Only tag path pushes five
   immutable GHCR tags and creates GitHub Release.
10. Download published assets; verify SHA256; run `--version`, `catalog`, `doctor`,
    example plotting, and one safe Recipe on clean machine. Record five GHCR
    digests and exact environment exports. Zenodo/JOSS are outside v0.4.0 scope.

Manual preflight with GitHub CLI:

```bash
gh workflow run release.yml --ref main -f release_tag=v0.4.0
gh run list --workflow release.yml --limit 5
gh run watch RUN_ID --exit-status
```

Annotated tag after successful preflight:

```bash
git fetch origin main
test "$(git rev-parse origin/main)" = "RECORDED_MAIN_SHA"
git tag -a v0.4.0 RECORDED_MAIN_SHA -m "BioHub v0.4.0"
git push origin v0.4.0
```
