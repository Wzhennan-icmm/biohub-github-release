#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


DNA = re.compile(r"^[ACGTNacgtn]+$")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait-id", required=True, action="append")
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--tail", required=True, choices=["5", "10"])
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()
    statuses = []
    candidates = []
    sequences = []
    failures = []
    for trait_id in args.trait_id:
        path = args.root / trait_id / "kmers" / f"pass_threshold_{args.tail}per"
        if not path.is_file():
            failures.append(f"{trait_id}\tmissing {path}")
            continue
        count = 0
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                count += 1
                fields = line.rstrip("\n").split()
                candidates.append((trait_id, line_number, line.rstrip("\n")))
                sequence = next((field.upper() for field in fields if DNA.fullmatch(field)), None)
                if sequence:
                    sequences.append((f"{trait_id}_kmer_{line_number}", sequence))
        statuses.append({"trait_id": trait_id, "significant_kmers": count, "source": str(path)})
    for path in [args.status, args.candidates, args.fasta, args.log]:
        path.parent.mkdir(parents=True, exist_ok=True)
    with args.status.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["trait_id", "significant_kmers", "source"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(statuses)
    with args.candidates.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait_id", "source_line", "raw_record"])
        writer.writerows(candidates)
    with args.fasta.open("w", encoding="utf-8") as handle:
        for identifier, sequence in sequences:
            handle.write(f">{identifier}\n{sequence}\n")
    args.log.write_text("\n".join(failures or ["all k-mer GWAS outputs parsed"]) + "\n", encoding="utf-8")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
