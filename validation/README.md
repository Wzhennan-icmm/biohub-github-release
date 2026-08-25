# BioHub validation register

`reviews.tsv` links each formerly or currently pending migration item to checked
evidence. `automated` means deterministic tests enforce evidence. `approved`
requires a named domain reviewer and review date. `pending` is never release-ready.

Synthetic tests do not replace biological validation. For domain, recipe, and
visual reviews, copy `review-template.md`, record deidentified dataset hashes,
reference software and versions, acceptance tolerances, reviewer, date, and
decision. Then update both migration matrix and `reviews.tsv` in same reviewed
commit.

Current human-review work is listed in
[`pending-review.zh-CN.md`](pending-review.zh-CN.md).

Four public CC0-1.0 representative packs live under `validation/packs/`:

- `annotation-coordinates`: inventory IDs 014, 015, 016
- `orthology-codon`: inventory IDs 030, 033, 034, 043
- `visualization`: inventory IDs 045, 046, 047
- `statistics`: inventory IDs 049, 050, 053, 072

Build self-contained evidence under ignored `validation/evidence/`:

```bash
python3 tools/validation_review.py build --pack all
python3 tools/validation_review.py verify --pack all
python3 tools/validation_review.py summary
```

`validation_review.py` has no approval operation. Generated `review.md` always
uses `human_status=pending`. Only explicit named human approval may change
`reviews.tsv` and `docs/SCRIPT_MIGRATION_MATRIX.tsv`.

Run normal consistency validation:

```bash
python3 tools/validate_release.py
```

Run formal release gate:

```bash
python3 tools/validate_release.py --release --tag v0.4.0
```
