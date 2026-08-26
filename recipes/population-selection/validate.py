#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def samples(path: Path):
    opener = gzip.open if path.suffix.lower() in {".gz", ".bgz"} else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                values = line.rstrip("\n").split("\t")[9:]
                if not values or any(not value for value in values) or len(values) != len(set(values)):
                    raise ValueError("VCF sample IDs must be non-empty and unique")
                return values
    raise ValueError("VCF has no sample header")


def sample_file(path: Path):
    values = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not values or len(values) != len(set(values)):
        raise ValueError(f"empty or duplicate sample IDs: {path}")
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    vcf_samples = set(samples(args.vcf))
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"comparison_id", "population1_samples", "population2_samples"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"manifest needs {sorted(required)}")
        rows = list(reader)
    normalized = []
    for row in rows:
        first_path = Path(row["population1_samples"]).expanduser().resolve()
        second_path = Path(row["population2_samples"]).expanduser().resolve()
        first = sample_file(first_path)
        second = sample_file(second_path)
        if set(first) & set(second):
            raise ValueError(f"overlapping populations for {row['comparison_id']}")
        missing = sorted((set(first) | set(second)) - vcf_samples)
        if missing:
            raise ValueError(f"VCF missing samples for {row['comparison_id']}: {missing[:10]}")
        normalized.append(
            {
                "comparison_id": row["comparison_id"],
                "population1_samples": str(first_path),
                "population2_samples": str(second_path),
                "population1_n": len(first),
                "population2_n": len(second),
            }
        )
    if not normalized or len({row["comparison_id"] for row in normalized}) != len(normalized):
        raise ValueError("comparison IDs must be non-empty and unique")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(normalized[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(normalized)
    args.summary.write_text(
        "metric\tvalue\n"
        f"vcf_samples\t{len(vcf_samples)}\n"
        f"comparisons\t{len(normalized)}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
