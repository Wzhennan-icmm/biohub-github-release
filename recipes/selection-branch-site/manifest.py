#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


ID_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")


def paml_header(path: Path) -> tuple[int, int]:
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if len(fields) < 2:
                break
            try:
                return int(fields[0]), int(fields[1])
            except ValueError:
                break
    raise ValueError(f"invalid PAML alignment header: {path}")


def load_manifest(path: Path, expected_taxa: int):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"test_id", "alignment", "marked_tree", "foreground"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"manifest must contain {sorted(required)}")
        rows = list(reader)
    if not rows:
        raise ValueError("tests manifest is empty")
    seen = set()
    normalized = []
    for row in rows:
        test_id = row["test_id"].strip()
        if not ID_PATTERN.fullmatch(test_id) or test_id in seen:
            raise ValueError(f"invalid or duplicate test_id: {test_id}")
        seen.add(test_id)
        alignment = Path(row["alignment"]).expanduser().resolve()
        tree = Path(row["marked_tree"]).expanduser().resolve()
        if not alignment.is_file() or not tree.is_file():
            raise FileNotFoundError(f"missing input for {test_id}: {alignment} or {tree}")
        tree_text = tree.read_text(encoding="utf-8").strip()
        if len(re.findall(r"#1(?![0-9])", tree_text)) != 1:
            raise ValueError(f"marked tree for {test_id} must contain exactly one #1")
        if not tree_text.endswith(";") or tree_text.count("(") != tree_text.count(")"):
            raise ValueError(f"marked tree for {test_id} is not a balanced semicolon-terminated Newick tree")
        foreground = row["foreground"].strip()
        if not foreground:
            raise ValueError(f"foreground label is empty for {test_id}")
        taxa, sites = paml_header(alignment)
        if taxa != expected_taxa:
            raise ValueError(f"{test_id}: alignment taxa {taxa} != expected {expected_taxa}")
        if sites <= 0 or sites % 3:
            raise ValueError(f"{test_id}: codon alignment length is not a positive multiple of three")
        normalized.append(
            {
                "test_id": test_id,
                "alignment": str(alignment),
                "marked_tree": str(tree),
                "foreground": foreground,
                "taxa": taxa,
                "sites": sites,
            }
        )
    return normalized


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--expected-taxa", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()
    rows = load_manifest(args.input, args.expected_taxa)
    for path in [args.output, args.summary, args.log]:
        path.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    args.summary.write_text(
        "metric\tvalue\n"
        f"tests\t{len(rows)}\n"
        f"expected_taxa\t{args.expected_taxa}\n"
        f"minimum_sites\t{min(row['sites'] for row in rows)}\n"
        f"maximum_sites\t{max(row['sites'] for row in rows)}\n",
        encoding="utf-8",
    )
    args.log.write_text("manifest validation passed\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
