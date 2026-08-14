#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def gamma_p(a: float, x: float) -> float:
    if x < 0 or a <= 0:
        raise ValueError("invalid incomplete-gamma arguments")
    if x == 0:
        return 0.0
    if x < a + 1:
        term = 1.0 / a
        total = term
        ap = a
        for _ in range(10000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-14:
                break
        return total * math.exp(-x + a * math.log(x) - math.lgamma(a))
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for index in range(1, 10000):
        an = -index * (index - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-14:
            break
    q = math.exp(-x + a * math.log(x) - math.lgamma(a)) * h
    return 1.0 - q


def chi_square_quantile(probability: float, degrees: float) -> float:
    if probability <= 0:
        return 0.0
    if probability >= 1:
        return math.inf
    low = 0.0
    high = max(1.0, degrees)
    while gamma_p(degrees / 2.0, high / 2.0) < probability:
        high *= 2.0
    for _ in range(200):
        middle = (low + high) / 2.0
        if gamma_p(degrees / 2.0, middle / 2.0) < probability:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def poisson_interval(count: int, confidence: float):
    alpha = 1.0 - confidence
    lower = 0.0 if count == 0 else 0.5 * chi_square_quantile(alpha / 2.0, 2.0 * count)
    upper = 0.5 * chi_square_quantile(1.0 - alpha / 2.0, 2.0 * (count + 1))
    return lower, upper


def read_bed(path: Path):
    intervals = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"invalid BED row at {path}:{line_number}")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"invalid BED interval at {path}:{line_number}")
            intervals.setdefault(chrom, []).append((start, end))
    merged = {}
    for chrom, values in intervals.items():
        result = []
        for start, end in sorted(values):
            if result and start <= result[-1][1]:
                result[-1] = (result[-1][0], max(result[-1][1], end))
            else:
                result.append((start, end))
        merged[chrom] = result
    if not merged:
        raise ValueError(f"callable BED is empty: {path}")
    return merged


def is_callable(intervals, chrom: str, position_one_based: int) -> bool:
    coordinate = position_one_based - 1
    for start, end in intervals.get(chrom, []):
        if coordinate < start:
            return False
        if start <= coordinate < end:
            return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--tier", required=True, action="append")
    parser.add_argument("--evidence-class", required=True, action="append")
    parser.add_argument("--ploidy", required=True, type=int)
    parser.add_argument("--confidence", required=True, type=float)
    parser.add_argument("--require-callable", required=True, choices=["true", "false"])
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--pair-rates", required=True, type=Path)
    parser.add_argument("--combined", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"family_id", "pair_id", "candidates_tsv", "callable_bed"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"pairs manifest needs {sorted(required)}")
        pairs = list(reader)
    if not pairs or len({row["pair_id"] for row in pairs}) != len(pairs):
        raise ValueError("pair_id values must be non-empty and unique")

    audit = []
    rates = []
    violations = []
    for pair in pairs:
        candidate_path = Path(pair["candidates_tsv"]).expanduser().resolve()
        bed_path = Path(pair["callable_bed"]).expanduser().resolve()
        intervals = read_bed(bed_path)
        callable_bp = sum(end - start for values in intervals.values() for start, end in values)
        with candidate_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"chrom", "position", "ref", "alt", "tier", "evidence_class"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(f"candidate table needs {sorted(required)}: {candidate_path}")
            candidates = list(reader)
        seen = set()
        included = 0
        for candidate in candidates:
            position = int(candidate["position"])
            key = (candidate["chrom"], position, candidate["ref"], candidate["alt"])
            if position < 1 or key in seen:
                raise ValueError(f"invalid or duplicate candidate in {candidate_path}: {key}")
            seen.add(key)
            callable_status = is_callable(intervals, candidate["chrom"], position)
            tier_status = candidate["tier"] in args.tier
            class_status = candidate["evidence_class"] in args.evidence_class
            use = callable_status and tier_status and class_status
            if not callable_status:
                violations.append(f"{pair['pair_id']}\t{candidate['chrom']}:{position}\tnot_callable")
            included += use
            audit.append(
                {
                    "family_id": pair["family_id"],
                    "pair_id": pair["pair_id"],
                    **candidate,
                    "callable": callable_status,
                    "tier_included": tier_status,
                    "evidence_class_included": class_status,
                    "rate_numerator_included": use,
                }
            )
        opportunity = callable_bp * args.ploidy
        lower_count, upper_count = poisson_interval(included, args.confidence)
        rates.append(
            {
                "family_id": pair["family_id"],
                "pair_id": pair["pair_id"],
                "candidate_count": included,
                "callable_bp": callable_bp,
                "ploidy": args.ploidy,
                "callable_opportunity": opportunity,
                "rate_per_bp_per_generation": included / opportunity,
                "ci_lower": lower_count / opportunity,
                "ci_upper": upper_count / opportunity,
                "confidence_level": args.confidence,
            }
        )
    if args.require_callable == "true" and violations:
        failure = True
    else:
        failure = False

    for path in [args.audit, args.pair_rates, args.combined, args.log]:
        path.parent.mkdir(parents=True, exist_ok=True)
    audit_fields = [
        "family_id", "pair_id", "chrom", "position", "ref", "alt", "tier", "evidence_class",
        "callable", "tier_included", "evidence_class_included", "rate_numerator_included",
    ]
    with args.audit.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(audit)
    rate_fields = list(rates[0])
    with args.pair_rates.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rate_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rates)
    total_count = sum(row["candidate_count"] for row in rates)
    total_opportunity = sum(row["callable_opportunity"] for row in rates)
    lower_count, upper_count = poisson_interval(total_count, args.confidence)
    args.combined.write_text(
        "candidate_count\tcallable_opportunity\trate_per_bp_per_generation\tci_lower\tci_upper\tconfidence_level\n"
        f"{total_count}\t{total_opportunity}\t{total_count / total_opportunity}\t"
        f"{lower_count / total_opportunity}\t{upper_count / total_opportunity}\t{args.confidence}\n",
        encoding="utf-8",
    )
    args.log.write_text("\n".join(violations or ["all included candidates are callable"]) + "\n", encoding="utf-8")
    return 1 if failure else 0


if __name__ == "__main__":
    raise SystemExit(main())
