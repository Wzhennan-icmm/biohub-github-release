#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


LNL = re.compile(r"lnL\([^\n]*?\):\s*([-+0-9.eE]+)")
BEB_LINE = re.compile(r"^\s*(\d+)\s+([A-Za-z*?-])\s+([0-9.]+)(\*{0,2})")


def parse_lnl(path: Path) -> float:
    text = path.read_text(encoding="utf-8", errors="replace")
    values = [float(match.group(1)) for match in LNL.finditer(text)]
    if not values:
        raise ValueError(f"cannot parse lnL from {path}")
    return values[-1]


def parse_beb(path: Path):
    text = path.read_text(encoding="utf-8", errors="replace")
    marker = "Bayes Empirical Bayes (BEB) analysis"
    if marker not in text:
        return []
    section = text.split(marker, 1)[1]
    records = []
    for line in section.splitlines():
        match = BEB_LINE.match(line)
        if match:
            records.append((int(match.group(1)), match.group(2), float(match.group(3)), match.group(4)))
    return records


def bh_adjust(values: list[float]) -> list[float]:
    count = len(values)
    order = sorted(range(count), key=values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for rank_index in range(count - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[original_index] * count / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--runs", required=True, type=Path)
    parser.add_argument("--beb-threshold", required=True, action="append", type=float)
    parser.add_argument("--lrt", required=True, type=Path)
    parser.add_argument("--beb", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        tests = list(csv.DictReader(handle, delimiter="\t"))
    results = []
    beb = []
    errors = []
    for test in tests:
        test_id = test["test_id"]
        alternative = args.runs / test_id / "alternative" / "mlc"
        null = args.runs / test_id / "null" / "mlc"
        try:
            alt_lnl = parse_lnl(alternative)
            null_lnl = parse_lnl(null)
            lrt = max(0.0, 2.0 * (alt_lnl - null_lnl))
            raw_p = 1.0 if lrt == 0 else 0.5 * math.erfc(math.sqrt(lrt / 2.0))
            results.append({**test, "alternative_lnL": alt_lnl, "null_lnL": null_lnl, "lrt": lrt, "pvalue_mixture": raw_p})
            for position, amino_acid, posterior, stars in parse_beb(alternative):
                for threshold in sorted(args.beb_threshold):
                    if posterior >= threshold:
                        beb.append(
                            {
                                "test_id": test_id,
                                "foreground": test["foreground"],
                                "position": position,
                                "amino_acid": amino_acid,
                                "posterior": posterior,
                                "threshold": threshold,
                                "paml_stars": stars or ".",
                            }
                        )
        except Exception as error:  # preserve remaining test summaries
            errors.append(f"{test_id}\t{error}")
    if results:
        adjusted = bh_adjust([row["pvalue_mixture"] for row in results])
        for row, qvalue in zip(results, adjusted):
            row["bh_qvalue_global"] = qvalue

    for path in [args.lrt, args.beb, args.summary, args.log]:
        path.parent.mkdir(parents=True, exist_ok=True)
    lrt_fields = [
        "test_id", "foreground", "alignment", "marked_tree", "taxa", "sites",
        "alternative_lnL", "null_lnL", "lrt", "pvalue_mixture", "bh_qvalue_global",
    ]
    with args.lrt.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=lrt_fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    beb_fields = ["test_id", "foreground", "position", "amino_acid", "posterior", "threshold", "paml_stars"]
    with args.beb.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=beb_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(beb)
    args.summary.write_text(
        "metric\tvalue\n"
        f"tests_requested\t{len(tests)}\n"
        f"tests_completed\t{len(results)}\n"
        f"tests_failed\t{len(errors)}\n"
        f"beb_records\t{len(beb)}\n",
        encoding="utf-8",
    )
    args.log.write_text("\n".join(errors or ["all model pairs parsed"]) + "\n", encoding="utf-8")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
