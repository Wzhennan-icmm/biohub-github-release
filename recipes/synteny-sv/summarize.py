#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-id", required=True, action="append")
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    failures = []
    for pair_id in args.pair_id:
        path = args.root / pair_id / "pair.syri.out"
        counts = Counter()
        records = 0
        if not path.is_file():
            failures.append(f"{pair_id}\tmissing {path}")
            continue
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    failures.append(
                        f"{pair_id}\tinvalid SyRI row {line_number}: expected at least 11 columns"
                    )
                    continue
                variant_type = fields[10].strip()
                if not variant_type:
                    failures.append(f"{pair_id}\tempty SyRI annotation type at row {line_number}")
                    continue
                counts[variant_type] += 1
                records += 1
        rows.append(
            {
                "pair_id": pair_id,
                "syri_records": records,
                "syri_type_counts": ";".join(f"{key}:{value}" for key, value in sorted(counts.items())),
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pair_id", "syri_records", "syri_type_counts"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    args.log.write_text("\n".join(failures or ["all SyRI pairs parsed"]) + "\n", encoding="utf-8")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
