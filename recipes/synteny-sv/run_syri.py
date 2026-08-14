#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--paf", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--query", required=True, type=Path)
    parser.add_argument("--cores", required=True, type=int)
    parser.add_argument("--include-cigar", action="store_true")
    parser.add_argument("--include-snps", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    command = [
        args.executable,
        "-c", str(args.paf.resolve()),
        "-r", str(args.reference.resolve()),
        "-q", str(args.query.resolve()),
        "-F", "P",
        "--prefix", "pair.",
        "--nc", str(args.cores),
        "--log", "INFO",
    ]
    if args.include_cigar:
        command.append("--cigar")
    if not args.include_snps:
        command.append("--nosnp")
    args.log.parent.mkdir(parents=True, exist_ok=True)
    with args.log.open("w", encoding="utf-8") as handle:
        handle.write("command=" + " ".join(command) + "\n")
        handle.flush()
        completed = subprocess.run(command, cwd=args.output, stdout=handle, stderr=subprocess.STDOUT, check=False)
    expected = args.output / "pair.syri.out"
    if completed.returncode == 0 and (not expected.is_file() or expected.stat().st_size == 0):
        raise RuntimeError("SyRI returned success without non-empty pair.syri.out")
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
