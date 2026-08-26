#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


def summarize_pi(path: Path, comparison_id: str, population: str):
    if not path.is_file():
        raise FileNotFoundError(path)
    values = []
    total_rows = 0
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "PI" not in reader.fieldnames:
            raise ValueError(f"windowed pi output lacks PI column: {path}")
        for row in reader:
            total_rows += 1
            raw = row.get("PI")
            if raw in {None, "nan", "-nan", "NA", ".", ""}:
                continue
            value = float(raw)
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid nucleotide diversity value in {path}: {raw}")
            values.append(value)
    return {
        "comparison_id": comparison_id,
        "population": population,
        "window_rows": total_rows,
        "finite_windows": len(values),
        "mean_window_pi": statistics.fmean(values) if values else "NA",
        "minimum_window_pi": min(values) if values else "NA",
        "maximum_window_pi": max(values) if values else "NA",
        "source_file": str(path),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison-id", required=True, action="append")
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--threshold", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--pi-summary", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()
    candidates = []
    summaries = []
    pi_summaries = []
    failures = []
    for comparison_id in args.comparison_id:
        path = args.root / comparison_id / "fst.windowed.weir.fst"
        if not path.is_file():
            failures.append(f"{comparison_id}\tmissing {path}")
            continue
        tested = 0
        selected = 0
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                value = row.get("MEAN_FST")
                if value in {None, "nan", "-nan", "NA", ""}:
                    continue
                fst = float(value)
                if not math.isfinite(fst):
                    continue
                tested += 1
                if fst >= args.threshold:
                    selected += 1
                    candidates.append(
                        {
                            "comparison_id": comparison_id,
                            "chromosome": row.get("CHROM", ""),
                            "window_start": row.get("BIN_START", ""),
                            "window_end": row.get("BIN_END", ""),
                            "variant_count": row.get("N_VARIANTS", ""),
                            "mean_fst": fst,
                        }
                    )
        summaries.append({"comparison_id": comparison_id, "tested_windows": tested, "candidate_windows": selected})
        for population in ["population1", "population2"]:
            try:
                pi_summaries.append(
                    summarize_pi(
                        args.root / comparison_id / f"{population}.windowed.pi",
                        comparison_id,
                        population,
                    )
                )
            except Exception as error:
                failures.append(f"{comparison_id}\t{population}\t{error}")
    for path in [args.output, args.summary, args.pi_summary, args.log]:
        path.parent.mkdir(parents=True, exist_ok=True)
    candidate_fields = ["comparison_id", "chromosome", "window_start", "window_end", "variant_count", "mean_fst"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=candidate_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(candidates)
    with args.summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["comparison_id", "tested_windows", "candidate_windows"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summaries)
    pi_fields = [
        "comparison_id", "population", "window_rows", "finite_windows",
        "mean_window_pi", "minimum_window_pi", "maximum_window_pi", "source_file",
    ]
    with args.pi_summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=pi_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(pi_summaries)
    args.log.write_text("\n".join(failures or ["all population comparisons parsed"]) + "\n", encoding="utf-8")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
