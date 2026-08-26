#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--genotype-kind", required=True, choices=["vcf", "pfile", "bfile"])
    parser.add_argument("--genotype", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--phenotype-column", required=True)
    parser.add_argument("--covariates")
    parser.add_argument("--covariate-columns")
    parser.add_argument("--maf", required=True, type=float)
    parser.add_argument("--geno", required=True, type=float)
    parser.add_argument("--mind", required=True, type=float)
    parser.add_argument("--hwe", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()

    genotype_flag = {"vcf": "--vcf", "pfile": "--pfile", "bfile": "--bfile"}[args.genotype_kind]
    command = [
        args.executable,
        genotype_flag,
        args.genotype,
        "--pheno", args.phenotype,
        "--pheno-name", args.phenotype_column,
        "--maf", str(args.maf),
        "--geno", str(args.geno),
        "--mind", str(args.mind),
        "--hwe", str(args.hwe),
        "--glm", "hide-covar", "allow-no-covars",
        "--out", str(args.output),
    ]
    if args.covariates:
        command.extend(["--covar", args.covariates])
        if args.covariate_columns:
            command.extend(["--covar-name", args.covariate_columns])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    with args.log.open("w", encoding="utf-8") as handle:
        handle.write("command=" + " ".join(command) + "\n")
        handle.flush()
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
