#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def fasta_records(path: Path):
    name = None
    chunks = []
    with path.open(encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"empty FASTA ID at line {line_number}")
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"sequence before header at line {line_number}")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def n50(lengths: list[int]) -> int:
    target = sum(lengths) / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= target:
            return length
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly-id", required=True)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--minimum-contig-length", required=True, type=int)
    parser.add_argument("--telomere-motif", required=True)
    parser.add_argument("--telomere-window", required=True, type=int)
    parser.add_argument("--minimum-telomere-hits", required=True, type=int)
    parser.add_argument("--sequences", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    motif = args.telomere_motif.upper()
    reverse = motif.translate(COMPLEMENT)[::-1]
    records = list(fasta_records(args.fasta))
    if not records or len({name for name, _ in records}) != len(records):
        raise ValueError("assembly FASTA is empty or contains duplicate IDs")
    rows = []
    for name, sequence in records:
        left = sequence[: args.telomere_window]
        right = sequence[-args.telomere_window :]
        left_hits = left.count(motif) + left.count(reverse)
        right_hits = right.count(motif) + right.count(reverse)
        rows.append(
            {
                "assembly_id": args.assembly_id,
                "sequence_id": name,
                "length": len(sequence),
                "n_bases": sequence.count("N"),
                "passes_minimum_length": len(sequence) >= args.minimum_contig_length,
                "left_telomere_hits": left_hits,
                "right_telomere_hits": right_hits,
                "both_ends_telomeric": left_hits >= args.minimum_telomere_hits
                and right_hits >= args.minimum_telomere_hits,
            }
        )
    args.sequences.parent.mkdir(parents=True, exist_ok=True)
    with args.sequences.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    lengths = [row["length"] for row in rows]
    args.summary.write_text(
        "assembly_id\tsequence_count\ttotal_bp\tN50\tlong_sequences\tboth_ends_telomeric\tn_bases\n"
        f"{args.assembly_id}\t{len(rows)}\t{sum(lengths)}\t{n50(lengths)}\t"
        f"{sum(row['passes_minimum_length'] for row in rows)}\t"
        f"{sum(row['both_ends_telomeric'] for row in rows)}\t"
        f"{sum(row['n_bases'] for row in rows)}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
