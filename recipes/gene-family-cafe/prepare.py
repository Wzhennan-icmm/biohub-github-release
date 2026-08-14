#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class Node:
    name: str = ""
    length: Optional[float] = None
    children: list["Node"] = field(default_factory=list)


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.index = 0

    def parse(self) -> Node:
        node = self._subtree()
        self._space()
        if self.index < len(self.text) and self.text[self.index] == ";":
            self.index += 1
        self._space()
        if self.index != len(self.text):
            raise ValueError(f"unexpected Newick text at position {self.index}")
        return node

    def _space(self) -> None:
        while self.index < len(self.text) and self.text[self.index].isspace():
            self.index += 1

    def _token(self) -> str:
        self._space()
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in ":,();":
            self.index += 1
        return self.text[start : self.index].strip()

    def _length(self) -> Optional[float]:
        self._space()
        if self.index >= len(self.text) or self.text[self.index] != ":":
            return None
        self.index += 1
        token = self._token()
        if not token:
            raise ValueError("empty branch length")
        value = float(token)
        if value <= 0:
            raise ValueError(f"non-positive branch length: {value}")
        return value

    def _subtree(self) -> Node:
        self._space()
        children: list[Node] = []
        if self.index < len(self.text) and self.text[self.index] == "(":
            self.index += 1
            while True:
                children.append(self._subtree())
                self._space()
                if self.index >= len(self.text):
                    raise ValueError("unterminated Newick subtree")
                delimiter = self.text[self.index]
                self.index += 1
                if delimiter == ")":
                    break
                if delimiter != ",":
                    raise ValueError(f"unexpected Newick delimiter: {delimiter}")
        name = self._token()
        length = self._length()
        return Node(name=name, length=length, children=children)


def leaf_depths(node: Node, depth: float = 0.0) -> dict[str, float]:
    here = depth + (node.length or 0.0)
    if not node.children:
        if not node.name:
            raise ValueError("unnamed leaf in species tree")
        return {node.name: here}
    if len(node.children) != 2:
        raise ValueError(f"species tree is not binary at node {node.name or '<unnamed>'}")
    result: dict[str, float] = {}
    for child in node.children:
        for name, value in leaf_depths(child, here).items():
            if name in result:
                raise ValueError(f"duplicate tree leaf: {name}")
            result[name] = value
    return result


def parse_counts(path: Path, input_format: str):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError("empty gene-count matrix") from error
        if input_format == "orthofinder":
            if not header or header[0] != "Orthogroup":
                raise ValueError("OrthoFinder matrix must start with Orthogroup")
            species_indices = [
                index for index, name in enumerate(header[1:], 1) if name != "Total"
            ]
            species = [header[index] for index in species_indices]
            id_index = 0
        else:
            normalized = [item.lstrip("#") for item in header]
            if "Family ID" in normalized:
                id_index = normalized.index("Family ID")
            elif "FamilyID" in normalized:
                id_index = normalized.index("FamilyID")
            else:
                raise ValueError("CAFE matrix needs Family ID or FamilyID column")
            species_indices = [
                index
                for index in range(id_index + 1, len(header))
                if header[index] != "Total"
            ]
            species = [header[index] for index in species_indices]
        if len(species) < 2 or len(set(species)) != len(species):
            raise ValueError("gene-count matrix needs at least two unique species columns")

        records = []
        seen = set()
        for line_number, row in enumerate(reader, 2):
            if not row or all(not cell.strip() for cell in row):
                continue
            if len(row) < len(header):
                raise ValueError(f"short row at line {line_number}")
            family = row[id_index].strip()
            if not family or family in seen:
                raise ValueError(f"empty or duplicate family ID at line {line_number}: {family}")
            seen.add(family)
            counts = []
            for index in species_indices:
                value = row[index].strip()
                if not value.isdigit():
                    raise ValueError(
                        f"non-negative integer required for {species[len(counts)]} at line {line_number}"
                    )
                counts.append(int(value))
            records.append((family, counts))
    if not records:
        raise ValueError("gene-count matrix has no families")
    return species, records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--input-format", required=True, choices=["orthofinder", "cafe"])
    parser.add_argument("--max-family-size", required=True, type=int)
    parser.add_argument("--ultrametric-tolerance", required=True, type=float)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--tree-output", required=True, type=Path)
    parser.add_argument("--filter-manifest", required=True, type=Path)
    parser.add_argument("--tree-qc", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    species, records = parse_counts(args.counts, args.input_format)
    tree_text = args.tree.read_text(encoding="utf-8").strip()
    depths = leaf_depths(NewickParser(tree_text).parse())
    if sorted(depths) != sorted(species):
        raise ValueError(
            f"tree/matrix taxa mismatch: tree={sorted(depths)} matrix={sorted(species)}"
        )
    deviation = max(depths.values()) - min(depths.values())
    if deviation > args.ultrametric_tolerance:
        raise ValueError(
            f"tree is not ultrametric: deviation {deviation} exceeds {args.ultrametric_tolerance}"
        )

    included = []
    manifest = []
    for family, counts in records:
        reasons = []
        if sum(counts) == 0:
            reasons.append("all_zero")
        if max(counts) > args.max_family_size:
            reasons.append("max_count_exceeds_limit")
        status = "excluded" if reasons else "included"
        manifest.append(
            [family, status, ";".join(reasons) or ".", str(max(counts)), str(sum(counts))]
        )
        if not reasons:
            included.append((family, counts))
    if not included:
        raise ValueError("all families were excluded")

    for path in [args.matrix, args.tree_output, args.filter_manifest, args.tree_qc, args.summary]:
        path.parent.mkdir(parents=True, exist_ok=True)
    with args.matrix.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Desc", "Family ID", *species])
        for family, counts in included:
            writer.writerow([family, family, *counts])
    args.tree_output.write_text(tree_text + "\n", encoding="utf-8")
    with args.filter_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["FamilyID", "status", "reason", "max_count", "total_count"])
        writer.writerows(manifest)
    with args.tree_qc.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxon", "root_to_tip"])
        writer.writerows(sorted(depths.items()))
        writer.writerow(["max_tip_deviation", deviation])
    args.summary.write_text(
        "metric\tvalue\n"
        f"input_families\t{len(records)}\n"
        f"included_families\t{len(included)}\n"
        f"excluded_families\t{len(records) - len(included)}\n"
        f"species\t{len(species)}\n"
        f"max_family_size\t{args.max_family_size}\n"
        f"max_tip_deviation\t{deviation}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
