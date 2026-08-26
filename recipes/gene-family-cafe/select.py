#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from pathlib import Path


PATTERNS = {
    "neg_log_likelihood": re.compile(r"Final Likelihood \(-lnL\):\s*([0-9.eE+-]+)"),
    "lambda": re.compile(r"Lambda:\s*([0-9.eE+-]+)"),
    "alpha": re.compile(r"Alpha:\s*([0-9.eE+-]+)"),
}


def parse_result(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    values = {}
    for name, pattern in PATTERNS.items():
        match = pattern.search(text)
        values[name] = float(match.group(1)) if match else math.nan
    values["numeric_warning"] = float("failure rates >20% of the time" in text)
    return values


def cv(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = abs(statistics.mean(values))
    return statistics.stdev(values) / mean if mean else math.inf


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", required=True, type=Path)
    parser.add_argument("--model-k", required=True, type=int)
    parser.add_argument("--replicates", required=True, type=int)
    parser.add_argument("--minimum-converged", required=True, type=int)
    parser.add_argument("--likelihood-tolerance", required=True, type=float)
    parser.add_argument("--parameter-cv-tolerance", required=True, type=float)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--selected", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()

    prefix = "Base" if args.model_k == 1 else "Gamma"
    records = []
    for replicate in range(1, args.replicates + 1):
        run_dir = args.runs / f"rep{replicate}"
        result = run_dir / f"{prefix}_results.txt"
        record = {"replicate": replicate, "run_dir": str(run_dir), "result_file": str(result)}
        if result.exists():
            record.update(parse_result(result))
        else:
            record.update(
                {"neg_log_likelihood": math.nan, "lambda": math.nan, "alpha": math.nan, "numeric_warning": 1.0}
            )
        records.append(record)
    valid = [
        row
        for row in records
        if math.isfinite(row["neg_log_likelihood"])
        and math.isfinite(row["lambda"])
        and not bool(row["numeric_warning"])
        and (args.model_k == 1 or math.isfinite(row["alpha"]))
    ]
    valid.sort(key=lambda row: row["neg_log_likelihood"])
    clusters = []
    for start in range(len(valid)):
        for stop in range(start + 1, len(valid) + 1):
            group = valid[start:stop]
            if group[-1]["neg_log_likelihood"] - group[0]["neg_log_likelihood"] > args.likelihood_tolerance:
                break
            if cv([row["lambda"] for row in group]) > args.parameter_cv_tolerance:
                continue
            if args.model_k > 1 and cv([row["alpha"] for row in group]) > args.parameter_cv_tolerance:
                continue
            clusters.append(group)
    if not clusters:
        selected_cluster = []
    else:
        selected_cluster = sorted(
            clusters,
            key=lambda group: (-len(group), statistics.median(row["neg_log_likelihood"] for row in group)),
        )[0]

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "replicate",
            "run_dir",
            "result_file",
            "neg_log_likelihood",
            "lambda",
            "alpha",
            "numeric_warning",
            "in_converged_cluster",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        selected_ids = {row["replicate"] for row in selected_cluster}
        for row in records:
            writer.writerow({**row, "in_converged_cluster": row["replicate"] in selected_ids})

    ok = len(selected_cluster) >= args.minimum_converged
    selected = min(selected_cluster, key=lambda row: row["neg_log_likelihood"]) if ok else None
    args.selected.write_text((selected["run_dir"] if selected else "") + "\n", encoding="utf-8")
    args.log.write_text(
        f"valid_replicates={len(valid)}\n"
        f"converged_replicates={len(selected_cluster)}\n"
        f"minimum_converged_replicates={args.minimum_converged}\n"
        f"selected_run={selected['run_dir'] if selected else ''}\n",
        encoding="utf-8",
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
