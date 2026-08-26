#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait-id", required=True, action="append")
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--threshold", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--lead", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()
    statuses = []
    leads = []
    failures = []
    for trait_id in args.trait_id:
        candidates = sorted((args.root / trait_id).glob("*.glm.*"))
        if not candidates:
            failures.append(f"{trait_id}\tno .glm result")
            continue
        tested = 0
        significant = 0
        best = None
        for path in candidates:
            with path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    normalized = {key.lstrip("#"): value for key, value in row.items()}
                    if normalized.get("TEST") not in {None, "ADD"}:
                        continue
                    raw_p = normalized.get("P")
                    if raw_p in {None, "NA", ".", ""}:
                        continue
                    pvalue = float(raw_p)
                    tested += 1
                    significant += pvalue <= args.threshold
                    if best is None or pvalue < best[0]:
                        best = (pvalue, normalized, path.name)
        statuses.append({"trait_id": trait_id, "result_files": len(candidates), "tested_rows": tested, "significant_rows": significant})
        if tested == 0:
            failures.append(f"{trait_id}\tno finite additive-test p-values")
        if best:
            pvalue, row, source = best
            leads.append(
                {
                    "trait_id": trait_id,
                    "chromosome": row.get("CHROM", ""),
                    "position": row.get("POS", ""),
                    "variant_id": row.get("ID", ""),
                    "pvalue": pvalue,
                    "passes_threshold": pvalue <= args.threshold,
                    "source_file": source,
                }
            )
    for path in [args.output, args.lead, args.log]:
        path.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["trait_id", "result_files", "tested_rows", "significant_rows"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(statuses)
    with args.lead.open("w", newline="", encoding="utf-8") as handle:
        fields = ["trait_id", "chromosome", "position", "variant_id", "pvalue", "passes_threshold", "source_file"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(leads)
    args.log.write_text("\n".join(failures or ["all trait results parsed"]) + "\n", encoding="utf-8")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
