#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation", required=True)
    parser.add_argument("--expected-taxa", required=True, type=int)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    source = Path(args.validation)
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"group", "protein_count", "cds_count", "status", "reason"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("validation_summary.tsv is empty or has invalid columns")

    completed = [row for row in rows if row["status"] == "completed"]
    failed = [row for row in rows if row["status"] != "completed"]
    mismatched = [
        row
        for row in completed
        if int(row["protein_count"]) != args.expected_taxa
        or int(row["cds_count"]) != args.expected_taxa
    ]

    summary = Path(args.summary)
    summary.parent.mkdir(parents=True, exist_ok=True)
    summary.write_text(
        "metric\tvalue\n"
        f"total_groups\t{len(rows)}\n"
        f"completed_groups\t{len(completed)}\n"
        f"failed_or_skipped_groups\t{len(failed)}\n"
        f"expected_taxa\t{args.expected_taxa}\n"
        f"taxon_count_mismatches\t{len(mismatched)}\n",
        encoding="utf-8",
    )

    log = Path(args.log)
    log.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"groups={len(rows)}",
        f"completed={len(completed)}",
        f"failed_or_skipped={len(failed)}",
        f"expected_taxa={args.expected_taxa}",
        f"taxon_count_mismatches={len(mismatched)}",
    ]
    for row in failed:
        lines.append(f"FAIL\t{row['group']}\t{row['status']}\t{row['reason']}")
    for row in mismatched:
        lines.append(
            "FAIL\t{}\tprotein_count={}\tcds_count={}".format(
                row["group"], row["protein_count"], row["cds_count"]
            )
        )
    log.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 1 if failed or mismatched else 0


if __name__ == "__main__":
    raise SystemExit(main())
