#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def fasta_ids(path: Path):
    identifiers = []
    seen = set()
    current = None
    has_sequence = False
    with path.open(encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, 1):
            value = line.strip()
            if not value:
                continue
            if value.startswith(">"):
                if current is not None and not has_sequence:
                    raise ValueError(f"empty FASTA record in {path}: {current}")
                identifier = value[1:].split()[0]
                if not identifier or identifier in seen:
                    raise ValueError(f"empty or duplicate FASTA ID at {path}:{line_number}")
                seen.add(identifier)
                identifiers.append(identifier)
                current = identifier
                has_sequence = False
            else:
                if current is None:
                    raise ValueError(f"sequence before FASTA header at {path}:{line_number}")
                has_sequence = True
    if current is not None and not has_sequence:
        raise ValueError(f"empty FASTA record in {path}: {current}")
    if not identifiers:
        raise ValueError(f"no FASTA records: {path}")
    return identifiers


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--query", required=True, type=Path)
    parser.add_argument("--require-matching-ids", required=True, choices=["true", "false"])
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    reference = fasta_ids(args.reference)
    query = fasta_ids(args.query)
    same = reference == query
    same_set = set(reference) == set(query)
    if args.require_matching_ids == "true" and not same_set:
        raise ValueError(
            f"SyRI pair has different sequence ID sets; reference-only={sorted(set(reference)-set(query))[:10]}, "
            f"query-only={sorted(set(query)-set(reference))[:10]}"
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "metric\tvalue\n"
        f"reference_sequences\t{len(reference)}\n"
        f"query_sequences\t{len(query)}\n"
        f"same_id_set\t{str(same_set).lower()}\n"
        f"same_id_order\t{str(same).lower()}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
