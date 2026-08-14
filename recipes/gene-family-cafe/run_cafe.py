#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--cores", required=True, type=int)
    parser.add_argument("--model-k", required=True, type=int)
    parser.add_argument("--root-distribution", required=True, type=float)
    parser.add_argument("--error-model")
    parser.add_argument("--lambda-value", type=float)
    parser.add_argument("--alpha-value", type=float)
    parser.add_argument("--family-pvalue", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    args = parser.parse_args()

    if args.model_k > 1 and args.alpha_value is None:
        raise SystemExit("--alpha-value is required when --model-k > 1")
    command = [
        args.executable,
        "-i",
        args.matrix,
        "-t",
        args.tree,
        "-c",
        str(args.cores),
        f"-p{args.root_distribution}",
        "-P",
        str(args.family_pvalue),
        "-o",
        str(args.output),
    ]
    if args.error_model is not None:
        command.append(f"-e{args.error_model}")
    if args.model_k > 1:
        command.extend(["-k", str(args.model_k), "-a", str(args.alpha_value)])
    if args.lambda_value is not None:
        command.extend(["-l", str(args.lambda_value)])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    with args.log.open("w", encoding="utf-8") as log:
        log.write("command=" + " ".join(command) + "\n")
        log.flush()
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
