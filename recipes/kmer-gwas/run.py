#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", required=True)
    parser.add_argument("--script", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--kmers-table", required=True)
    parser.add_argument("--kmers-number", required=True, type=int)
    parser.add_argument("--permutations", required=True, type=int)
    parser.add_argument("--maf", required=True, type=float)
    parser.add_argument("--mac", required=True, type=int)
    parser.add_argument("--minimum-data-points", required=True, type=int)
    parser.add_argument("--kmer-length", required=True, type=int)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        args.python,
        args.script,
        "--min_data_points", str(args.minimum_data_points),
        "--pheno", args.phenotype,
        "--kmers_table", args.kmers_table,
        "--kmers_number", str(args.kmers_number),
        "--permutations", str(args.permutations),
        "--maf", str(args.maf),
        "--mac", str(args.mac),
        "-l", str(args.kmer_length),
        "-p", str(args.threads),
        "--outdir", str(args.output),
    ]
    args.log.parent.mkdir(parents=True, exist_ok=True)
    with args.log.open("w", encoding="utf-8") as handle:
        handle.write("command=" + " ".join(command) + "\n")
        handle.flush()
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
