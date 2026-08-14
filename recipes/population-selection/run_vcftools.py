#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def input_args(vcf: str):
    return ["--gzvcf", vcf] if vcf.lower().endswith((".gz", ".bgz")) else ["--vcf", vcf]


def run(command, log: Path):
    with log.open("a", encoding="utf-8") as handle:
        handle.write("command=" + " ".join(command) + "\n")
        handle.flush()
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if completed.returncode:
        raise SystemExit(completed.returncode)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--pop1", required=True)
    parser.add_argument("--pop2", required=True)
    parser.add_argument("--window", required=True, type=int)
    parser.add_argument("--step", required=True, type=int)
    parser.add_argument("--maf", required=True, type=float)
    parser.add_argument("--max-missing", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    common = [
        args.executable,
        *input_args(args.vcf),
        "--maf", str(args.maf),
        "--max-missing", str(args.max_missing),
    ]
    run(
        [
            *common,
            "--weir-fst-pop", args.pop1,
            "--weir-fst-pop", args.pop2,
            "--fst-window-size", str(args.window),
            "--fst-window-step", str(args.step),
            "--out", str(args.output / "fst"),
        ],
        args.log,
    )
    for label, population in [("population1", args.pop1), ("population2", args.pop2)]:
        run(
            [
                *common,
                "--keep", population,
                "--window-pi", str(args.window),
                "--window-pi-step", str(args.step),
                "--out", str(args.output / label),
            ],
            args.log,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
