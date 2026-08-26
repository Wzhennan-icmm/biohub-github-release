# Script migration ledger

BioHub separates private discovery evidence from publishable migration decisions.
`SCRIPT_MIGRATION_MATRIX.tsv` contains sanitized source names and one decision for
every candidate script visible in workspace audit on 2026-08-14. It does not copy
unpublished data, absolute paths, credentials, sample identifiers, or script bodies.

## Rebuild private inventory

Run from repository root. Output path is ignored by Git:

```bash
python3 tools/audit_script_inventory.py \
  --root .. \
  --exclude biohub-github-release/recipes \
  --exclude biohub-github-release/compat \
  --exclude biohub-github-release/tools \
  --exclude biohub-github-release/r/psmc_plot.R \
  --output audit/private-script-inventory.tsv
```

Exclusions remove newly integrated recipe/tool targets from source-candidate
discovery; retained pre-migration implementations remain in the scan. If source
workspace gains new scripts, append matrix rows instead of silently adding another
exclusion.

Review private checksums against matrix before changing migration state. Never
commit private inventory: it contains absolute paths. A repeated scan is
idempotent when inputs are unchanged; replacement requires explicit `--force`.

## Decision meanings

- `integrated`: callable BioHub command has replaced script interface.
- `recipe`: parameterized, provenance-producing workflow covers scientific intent.
- `superseded`: duplicate launcher/module replaced by existing BioHub surface.
- `retained`: file is current maintained implementation or launcher.
- `deferred`: useful candidate, but inputs, license, method, or golden outputs remain unresolved.
- `excluded`: tutorial, malformed notebook-like code, interactive app, or third-party code unsuitable for copying into core.

`integrated` and `recipe` do not mean publication validation is complete. Matrix
verification column distinguishes syntax/CLI coverage from domain-approved golden
outputs. Advanced methods remain experimental until representative deidentified
fixtures and expert-approved expected results exist.

Release evidence is tracked in `../validation/reviews.tsv`. `automated` entries
are enforced by deterministic golden tests. `approved` requires a named human
reviewer and ISO review date. `pending` blocks formal tag publication through
`tools/validate_release.py --release`.

## 2026-08-14 snapshot

| Decision | Count |
| --- | ---: |
| integrated | 63 |
| recipe | 5 |
| superseded | 11 |
| retained | 2 |
| deferred | 5 |
| excluded | 5 |
| total | 91 |

Thus 68 source candidates now map to callable commands or recipes. Deferred and
excluded rows remain visible with reasons; they were not copied into core by
guessing missing methods, licenses, inputs, or expected outputs.

## Migration contract

1. Preserve source checksum and decision before edits.
2. Add new command/recipe without removing legacy entry.
3. Validate inputs, explicit parameters, output boundaries, and overwrite behavior.
4. Compare sanitized fixture outputs; record known differences.
5. Update catalog, manual, CI, and matrix in same change.
6. Remove compatibility path only in announced major-version window.

## Deferred v0.5 candidates

Deferred rows stay outside v0.4 release scope. Promotion requires all four entry
criteria: redistributable source/license, explicit input and output schema,
parameterized method without project-specific paths or runtime package installs,
and a sanitized golden fixture with named reviewer acceptance.

| Candidate | Current blocker | v0.5 entry evidence |
| --- | --- | --- |
| `Chinese_map.R` | External map services, shapefiles, fonts, non-genomics scope | Pinned redistributable map assets, offline render contract, visual golden |
| `Chromosome_plot.R` | Runtime installs and undefined project objects | Declared table schema, fixed dependencies, coordinate and visual goldens |
| `Pheatmap.R` | Hard-coded exploratory plot | Reusable heatmap schema, scaling/clustering parameters, visual golden |
| `gggenes.R` | Package-installing examples, no approved contract | Fixed dependencies, feature schema, layout parameters, visual golden |
| `hamstr.sh` | Hard-coded taxa/paths, unsafe iteration, license unresolved | License approval, manifest-driven inputs, safe runner, end-to-end golden |

Opening a v0.5 implementation issue must name owner, reviewer, fixture hashes,
acceptance tolerance, and target command or recipe. Missing evidence keeps row
`deferred`; it must not be copied into core speculatively.
