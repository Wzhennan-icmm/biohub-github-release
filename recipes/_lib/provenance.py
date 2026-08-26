#!/usr/bin/env python3
"""Small standard-library helpers shared by packaged BioHub recipes."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
from typing import Iterable, Sequence, Tuple


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _regular_files(path: Path) -> Iterable[Path]:
    if path.is_file() and not path.is_symlink():
        yield path
        return
    if path.is_dir():
        for child in sorted(path.rglob("*")):
            if child.is_file() and not child.is_symlink():
                yield child


def write_input_manifest(
    items: Sequence[Tuple[str, str]], output: str | os.PathLike[str]
) -> None:
    rows = ["logical_name\tpath\tsize_bytes\tsha256"]
    for logical_name, raw_path in items:
        path = Path(raw_path).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"missing input for {logical_name}: {path}")
        files = list(_regular_files(path))
        if not files:
            raise ValueError(f"input contains no regular files: {path}")
        for file_path in files:
            suffix = file_path.name if path.is_file() else str(file_path.relative_to(path))
            name = logical_name if path.is_file() else f"{logical_name}/{suffix}"
            rows.append(
                f"{name}\t{file_path}\t{file_path.stat().st_size}\t{sha256_file(file_path)}"
            )
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_readme(
    output: str | os.PathLike[str],
    recipe_id: str,
    purpose: str,
    validation: Sequence[str],
) -> None:
    lines = [
        f"BioHub recipe: {recipe_id}",
        "",
        purpose,
        "",
        "Validation:",
    ]
    lines.extend(f"- {item}" for item in validation)
    lines.extend(
        [
            "",
            "See config.resolved.yaml, command.sh, versions.tsv, provenance.json,",
            "inputs.manifest.tsv, logs/, recipe-declared result tables and checksums.sha256.",
        ]
    )
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
