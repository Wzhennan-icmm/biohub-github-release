#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from pathlib import Path


def autocorrelation(values: list[float], lag: int, mean: float, variance: float) -> float:
    count = len(values) - lag
    if count <= 0 or variance == 0:
        return 0.0
    covariance = sum((values[i] - mean) * (values[i + lag] - mean) for i in range(count)) / count
    return covariance / variance


def effective_sample_size(values: list[float]) -> float:
    count = len(values)
    if count < 3:
        return float(count)
    mean = statistics.mean(values)
    variance = statistics.pvariance(values)
    if variance == 0:
        return float(count)
    rho_sum = 0.0
    lag = 1
    while lag + 1 < count:
        pair = autocorrelation(values, lag, mean, variance) + autocorrelation(values, lag + 1, mean, variance)
        if pair <= 0:
            break
        rho_sum += pair
        lag += 2
    return count / (1.0 + 2.0 * rho_sum)


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def read_chain(path: Path, burnin: int, pattern: re.Pattern[str]):
    with path.open(encoding="utf-8") as handle:
        rows = [line.split() for line in handle if line.strip()]
    if len(rows) <= burnin + 1:
        raise ValueError(f"chain has too few samples after burn-in: {path}")
    header = rows[0]
    if len(header) != len(set(header)):
        raise ValueError(f"chain header has duplicate columns: {path}")
    indices = [index for index, name in enumerate(header) if pattern.search(name)]
    if not indices:
        raise ValueError(f"no age columns match configured regex in {path}")
    data = {header[index]: [] for index in indices}
    for row_number, row in enumerate(rows[burnin + 1 :], burnin + 2):
        if len(row) != len(header):
            raise ValueError(
                f"chain row {row_number} has {len(row)} fields; expected {len(header)}: {path}"
            )
        for index in indices:
            value = float(row[index])
            if not math.isfinite(value):
                raise ValueError(f"non-finite chain value at row {row_number}: {path}")
            data[header[index]].append(value)
    return data


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", required=True, type=Path)
    parser.add_argument("--run-id", required=True, action="append")
    parser.add_argument("--age-column-regex", required=True)
    parser.add_argument("--burnin-samples", required=True, type=int)
    parser.add_argument("--expected-nodes", required=True, type=int)
    parser.add_argument("--minimum-ess", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()

    if args.burnin_samples < 0 or args.expected_nodes < 1 or args.minimum_ess <= 0:
        raise ValueError("burn-in, expected nodes, or minimum ESS is outside valid range")

    pattern = re.compile(args.age_column_regex)
    records = []
    failures = []
    for run_id in args.run_id:
        try:
            chain = read_chain(args.runs / run_id / "mcmc.txt", args.burnin_samples, pattern)
            if len(chain) != args.expected_nodes:
                raise ValueError(f"matched {len(chain)} nodes, expected {args.expected_nodes}")
            for node, values in chain.items():
                ess = effective_sample_size(values)
                records.append(
                    {
                        "run_id": run_id,
                        "node": node,
                        "samples": len(values),
                        "mean": statistics.mean(values),
                        "sd": statistics.stdev(values) if len(values) > 1 else 0.0,
                        "median": statistics.median(values),
                        "minimum": min(values),
                        "maximum": max(values),
                        "ci_lower_2.5pct": quantile(values, 0.025),
                        "ci_upper_97.5pct": quantile(values, 0.975),
                        "ess": ess,
                        "ess_pass": ess >= args.minimum_ess,
                    }
                )
        except Exception as error:
            failures.append(f"{run_id}\t{error}")

    fields = [
        "run_id", "node", "samples", "mean", "sd", "median", "minimum", "maximum",
        "ci_lower_2.5pct", "ci_upper_97.5pct", "ess", "ess_pass",
    ]
    for path in [args.output, args.summary, args.log]:
        path.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    failed_ess = sum(not row["ess_pass"] for row in records)
    args.summary.write_text(
        "metric\tvalue\n"
        f"runs_requested\t{len(args.run_id)}\n"
        f"runs_failed\t{len(failures)}\n"
        f"node_run_records\t{len(records)}\n"
        f"records_below_minimum_ess\t{failed_ess}\n"
        f"minimum_ess\t{args.minimum_ess}\n",
        encoding="utf-8",
    )
    args.log.write_text("\n".join(failures or ["all chains parsed"]) + "\n", encoding="utf-8")
    return 1 if failures or failed_ess else 0


if __name__ == "__main__":
    raise SystemExit(main())
