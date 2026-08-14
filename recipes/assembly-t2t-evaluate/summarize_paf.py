#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def union_length(intervals):
    total = 0
    end = -1
    for start, stop in sorted(intervals):
        if stop <= end:
            continue
        total += stop - max(start, end)
        end = stop
    return total


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    identifier = None
    length = 0
    with path.open(encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, 1):
            value = line.strip()
            if not value:
                continue
            if value.startswith(">"):
                if identifier is not None:
                    if length == 0:
                        raise ValueError(f"empty FASTA record: {path}:{identifier}")
                    lengths[identifier] = length
                identifier = value[1:].split()[0]
                if not identifier or identifier in lengths:
                    raise ValueError(f"empty or duplicate FASTA ID at {path}:{line_number}")
                length = 0
            else:
                if identifier is None:
                    raise ValueError(f"sequence before FASTA header at {path}:{line_number}")
                length += len(value)
    if identifier is not None:
        if length == 0:
            raise ValueError(f"empty FASTA record: {path}:{identifier}")
        lengths[identifier] = length
    if not lengths:
        raise ValueError(f"FASTA contains no records: {path}")
    return lengths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly-id", required=True)
    parser.add_argument("--paf", required=True, type=Path)
    parser.add_argument("--query-fasta", required=True, type=Path)
    parser.add_argument("--target-fasta", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    query_intervals = defaultdict(list)
    target_intervals = defaultdict(list)
    query_lengths = fasta_lengths(args.query_fasta)
    target_lengths = fasta_lengths(args.target_fasta)
    alignments = 0
    matches = 0
    block_bases = 0
    with args.paf.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"invalid PAF row at line {line_number}")
            qname, qlen, qstart, qend = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            tname, tlen, tstart, tend = fields[5], int(fields[6]), int(fields[7]), int(fields[8])
            if qname not in query_lengths or tname not in target_lengths:
                raise ValueError(f"PAF row {line_number} names a sequence absent from FASTA")
            if qlen != query_lengths[qname] or tlen != target_lengths[tname]:
                raise ValueError(f"PAF/FASTA length mismatch at line {line_number}")
            if not (0 <= qstart < qend <= qlen and 0 <= tstart < tend <= tlen):
                raise ValueError(f"invalid PAF interval at line {line_number}")
            row_matches = int(fields[9])
            row_block_bases = int(fields[10])
            if row_matches < 0 or row_block_bases <= 0 or row_matches > row_block_bases:
                raise ValueError(f"invalid PAF match counts at line {line_number}")
            query_intervals[qname].append((qstart, qend))
            target_intervals[tname].append((tstart, tend))
            matches += row_matches
            block_bases += row_block_bases
            alignments += 1
    if alignments == 0:
        raise ValueError("PAF contains no alignments")
    query_aligned = sum(union_length(value) for value in query_intervals.values())
    target_aligned = sum(union_length(value) for value in target_intervals.values())
    query_total = sum(query_lengths.values())
    target_total = sum(target_lengths.values())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "assembly_id", "alignments", "query_sequences", "target_sequences",
                "query_union_coverage", "target_union_coverage", "weighted_identity",
            ]
        )
        writer.writerow(
            [
                args.assembly_id,
                alignments,
                len(query_lengths),
                len(target_lengths),
                query_aligned / query_total if query_total else 0,
                target_aligned / target_total if target_total else 0,
                matches / block_bases if block_bases else 0,
            ]
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
