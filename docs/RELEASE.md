# Release checklist

1. All CI jobs green on `main`; every stable command has approved fixtures.
2. Update version, `CHANGELOG.md`, `CITATION.cff`, docs, and compatibility notes.
3. Build Linux/macOS binaries; publish SHA256 checksums and Docker image digest.
4. Create annotated `vX.Y.Z` tag; GitHub Actions publishes release assets.
5. Archive tagged release with Zenodo and add generated DOI to `CITATION.cff`.
6. Record known limitations, external dependency versions, and tested platforms.
7. For JOSS: verify public-history timing, paper metadata, impact citations, and AI disclosure.
8. Run `python3 tools/validate_recipes.py --require-r`; verify every recipe path is packaged.
9. For domain images, record immutable image digest and exact environment export; never publish a mutable-only tag.
10. Run recipe helper regressions, all 13 DAG dry-runs, and four synthetic safe executions; archive logs.
11. Verify release archive contains `docs/USER_GUIDE.zh-CN.md` and `docs/RECIPES.zh-CN.md`; validate all local Markdown links.
