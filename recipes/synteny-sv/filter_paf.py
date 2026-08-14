#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--minimum-length", required=True, type=int)
    parser.add_argument("--minimum-mapq", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    kept = 0
    total = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.input.open(encoding="utf-8") as source, args.output.open("w", encoding="utf-8") as destination:
        for line_number, line in enumerate(source, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"invalid PAF row at line {line_number}")
            total += 1
            if int(fields[10]) >= args.minimum_length and int(fields[11]) >= args.minimum_mapq:
                destination.write(line)
                kept += 1
    if total == 0 or kept == 0:
        raise ValueError(f"PAF filtering retained {kept} of {total} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
