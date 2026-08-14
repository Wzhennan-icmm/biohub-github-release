#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def control_text(args, seqfile: Path, treefile: Path, outfile: Path) -> str:
    is_null = args.model == "null"
    omega = 1.0 if is_null else args.alternative_initial_omega
    return f"""seqfile = {seqfile}
treefile = {treefile}
outfile = {outfile}
noisy = 3
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = {args.codon_frequency}
clock = 0
aaDist = 0
model = 2
NSsites = 2
icode = 0
Mgene = 0
fix_kappa = 0
kappa = {args.initial_kappa}
fix_omega = {1 if is_null else 0}
omega = {omega}
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8
getSE = 0
RateAncestor = 0
Small_Diff = 1e-6
cleandata = {args.clean_data}
method = {args.optimizer_method}
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--model", required=True, choices=["alternative", "null"])
    parser.add_argument("--codon-frequency", required=True, type=int)
    parser.add_argument("--initial-kappa", required=True, type=float)
    parser.add_argument("--alternative-initial-omega", required=True, type=float)
    parser.add_argument("--clean-data", required=True, type=int)
    parser.add_argument("--optimizer-method", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=False)
    sequence = args.output / "alignment.paml"
    tree = args.output / "marked_tree.nwk"
    shutil.copyfile(args.alignment, sequence)
    shutil.copyfile(args.tree, tree)
    control = args.output / "codeml.ctl"
    mlc = args.output / "mlc"
    control.write_text(control_text(args, sequence.resolve(), tree.resolve(), mlc.resolve()), encoding="utf-8")
    args.log.parent.mkdir(parents=True, exist_ok=True)
    with args.log.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            [args.executable, str(control.resolve())],
            cwd=args.output,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
