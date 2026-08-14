#!/usr/bin/env python3
"""Create private, checksum-backed script inventory for migration audits."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile


LANGUAGES = {
    ".awk": "awk",
    ".pl": "perl",
    ".py": "python",
    ".r": "r",
    ".sh": "shell",
    ".smk": "snakemake",
}
DEFAULT_EXCLUDED_NAMES = {".git", ".venv", "node_modules", "target", "venv", "__pycache__"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inventory scripts without copying their contents. Output can contain private paths."
    )
    parser.add_argument("--root", required=True, help="Directory to scan")
    parser.add_argument("--output", required=True, help="Private .tsv or .json output")
    parser.add_argument("--exclude", action="append", default=[], help="Path under root to exclude; repeatable")
    parser.add_argument("--force", action="store_true", help="Replace existing output")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def is_excluded(path: Path, root: Path, excluded: list[Path]) -> bool:
    relative = path.relative_to(root)
    if any(part in DEFAULT_EXCLUDED_NAMES for part in relative.parts):
        return True
    return any(path == item or item in path.parents for item in excluded)


def collect(root: Path, excluded: list[Path]) -> list[dict[str, object]]:
    rows = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.is_symlink() or is_excluded(path, root, excluded):
            continue
        suffix = path.suffix.lower()
        language = "snakemake" if path.name == "Snakefile" else LANGUAGES.get(suffix)
        if language is None:
            continue
        rows.append(
            {
                "relative_path": str(path.relative_to(root)),
                "absolute_path": str(path),
                "language": language,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return rows


def serialize(rows: list[dict[str, object]], suffix: str) -> str:
    if suffix == ".json":
        return json.dumps({"schema_version": 1, "files": rows}, ensure_ascii=False, indent=2) + "\n"
    header = ["relative_path", "absolute_path", "language", "size_bytes", "sha256"]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=header, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def main() -> int:
    args = parse_args()
    root = Path(args.root).expanduser().resolve()
    output = Path(args.output).expanduser().resolve()
    if not root.is_dir():
        raise SystemExit(f"root is not a directory: {root}")
    if output.exists() and not args.force:
        raise SystemExit(f"output exists: {output}; pass --force to replace")
    if output.suffix.lower() not in {".tsv", ".json"}:
        raise SystemExit("output extension must be .tsv or .json")
    excluded = [(root / value).resolve() for value in args.exclude]
    if any(root != item and root not in item.parents for item in excluded):
        raise SystemExit("excluded paths must stay under root")
    rows = collect(root, excluded)
    output.parent.mkdir(parents=True, exist_ok=True)
    content = serialize(rows, output.suffix.lower())
    handle = tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=output.parent, delete=False)
    temporary = Path(handle.name)
    try:
        with handle:
            handle.write(content)
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()
    print(f"inventoried={len(rows)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
