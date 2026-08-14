#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_one(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"expected one summary row: {path}")
    return rows[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly-id", required=True, action="append")
    parser.add_argument("--stats-dir", required=True, type=Path)
    parser.add_argument("--paf-dir", required=True, type=Path)
    parser.add_argument("--expected-chromosomes", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    warnings = []
    for assembly_id in args.assembly_id:
        stats = read_one(args.stats_dir / f"{assembly_id}.summary.tsv")
        mapping = read_one(args.paf_dir / f"{assembly_id}.tsv")
        row = {**stats, **{key: value for key, value in mapping.items() if key != "assembly_id"}}
        row["expected_chromosomes"] = args.expected_chromosomes
        row["sequence_count_matches_expected"] = int(stats["sequence_count"]) == args.expected_chromosomes
        if not row["sequence_count_matches_expected"]:
            warnings.append(
                f"{assembly_id}: sequence_count={stats['sequence_count']} expected={args.expected_chromosomes}"
            )
        rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    args.log.write_text("\n".join(warnings or ["all assembly summaries passed configured checks"]) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
