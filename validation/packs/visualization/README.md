# Visualization validation pack

Public CC0-1.0 depth fixtures cover inventory IDs `045`, `046`, and `047`.
Automated checks lock transformed tables, titles, point counts, finite SVG
coordinates, and whole-file FNV-1a fingerprints.

Build and verify:

```bash
python3 tools/validation_review.py build --pack visualization
python3 tools/validation_review.py verify --pack visualization
open validation/evidence/visualization/gallery.html
```

Fingerprint equality does not replace visual review. Open `gallery.html`; inspect
at 100% and at intended publication-column size before approval.
