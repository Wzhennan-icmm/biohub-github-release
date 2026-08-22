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

Run normal consistency validation:

```bash
python3 tools/validate_release.py
```

Run formal release gate:

```bash
python3 tools/validate_release.py --release --tag v0.4.0
```
