#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected", required=True, type=Path)
    parser.add_argument("--model-k", required=True, type=int)
    parser.add_argument("--expanded", required=True, type=Path)
    parser.add_argument("--contracted", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    run_dir = Path(args.selected.read_text(encoding="utf-8").strip())
    if not run_dir:
        raise SystemExit("selected CAFE run is empty")
    prefix = "Base" if args.model_k == 1 else "Gamma"
    path = run_dir / f"{prefix}_change.tab"
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if len(rows) < 2 or len(rows[0]) < 2:
        raise ValueError("invalid CAFE change table")
    header = [value.lstrip("#") for value in rows[0]]
    family_index = header.index("FamilyID") if "FamilyID" in header else 0
    changes = []
    counts: Counter[tuple[str, str]] = Counter()
    for row in rows[1:]:
        if not row:
            continue
        family = row[family_index]
        for index, node in enumerate(header):
            if index == family_index or index >= len(row) or not row[index].strip():
                continue
            value = int(float(row[index]))
            if value == 0:
                continue
            direction = "expansion" if value > 0 else "contraction"
            changes.append((direction, family, node, value))
            counts[(node, direction)] += 1
    for destination, direction in [(args.expanded, "expansion"), (args.contracted, "contraction")]:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["FamilyID", "node", "change"])
            writer.writerows((family, node, value) for kind, family, node, value in changes if kind == direction)
    with args.summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["node", "direction", "families"])
        for (node, direction), count in sorted(counts.items()):
            writer.writerow([node, direction, count])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
