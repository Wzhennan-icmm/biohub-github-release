#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path


SETTING = re.compile(r"^(\s*)([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*?)(\s*(?:\*.*)?)$")


def parse_settings(path: Path) -> dict[str, str]:
    settings = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = SETTING.match(line)
        if match:
            settings[match.group(2).lower()] = match.group(3).strip()
    return settings


def first_token(value: str) -> str:
    return value.split()[0] if value.split() else ""


def resolve_input(control: Path, value: str) -> Path:
    path = Path(first_token(value)).expanduser()
    if not path.is_absolute():
        path = control.parent / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def rewrite_control(source: Path, destination: Path, output: Path, in_bv: Path | None) -> None:
    settings = parse_settings(source)
    replacements = {
        "seqfile": str(resolve_input(source, settings["seqfile"])),
        "treefile": str(resolve_input(source, settings["treefile"])),
        "outfile": str(output.resolve()),
    }
    if in_bv is not None:
        tokens = settings["usedata"].split()
        replacements["usedata"] = " ".join([tokens[0], str(in_bv.resolve()), *tokens[2:]])
    lines = []
    for line in source.read_text(encoding="utf-8", errors="replace").splitlines():
        match = SETTING.match(line)
        if match and match.group(2).lower() in replacements:
            key = match.group(2)
            lines.append(f"{key} = {replacements[key.lower()]}")
        else:
            lines.append(line)
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_controls(stage1: Path, stage2: Path, expected_loci: int) -> None:
    first = parse_settings(stage1)
    second = parse_settings(stage2)
    required = {"seqfile", "treefile", "outfile", "ndata", "usedata"}
    for name, settings in [("stage1", first), ("stage2", second)]:
        missing = required - settings.keys()
        if missing:
            raise ValueError(f"{name} control missing settings: {sorted(missing)}")
        if int(first_token(settings["ndata"])) != expected_loci:
            raise ValueError(f"{name} ndata does not equal expected_loci={expected_loci}")
        resolve_input(stage1 if name == "stage1" else stage2, settings["seqfile"])
        resolve_input(stage1 if name == "stage1" else stage2, settings["treefile"])
    if first_token(first["usedata"]) != "3":
        raise ValueError("stage1 usedata must equal 3")
    if first_token(second["usedata"]) != "2":
        raise ValueError("stage2 usedata must equal 2")


def run(executable: str, control: Path, cwd: Path, log: Path) -> int:
    with log.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(
            [executable, str(control.resolve())],
            cwd=cwd,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return completed.returncode


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--stage1", required=True, type=Path)
    parser.add_argument("--stage2", required=True, type=Path)
    parser.add_argument("--expected-loci", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log-prefix", required=True, type=Path)
    args = parser.parse_args()

    stage1 = args.stage1.expanduser().resolve()
    stage2 = args.stage2.expanduser().resolve()
    validate_controls(stage1, stage2, args.expected_loci)
    args.output.mkdir(parents=True, exist_ok=False)
    args.log_prefix.parent.mkdir(parents=True, exist_ok=True)

    first_ctl = args.output / "stage1.usedata3.ctl"
    second_ctl = args.output / "stage2.usedata2.ctl"
    rewrite_control(stage1, first_ctl, args.output / "stage1.mcmctree.out", None)
    if run(args.executable, first_ctl, args.output, Path(str(args.log_prefix) + ".stage1.log")):
        return 1
    out_bv = args.output / "out.BV"
    if not out_bv.is_file() or out_bv.stat().st_size == 0:
        raise RuntimeError("stage1 did not produce non-empty out.BV")
    in_bv = args.output / "in.BV"
    shutil.copyfile(out_bv, in_bv)
    rewrite_control(stage2, second_ctl, args.output / "mcmctree.out", in_bv)
    if run(args.executable, second_ctl, args.output, Path(str(args.log_prefix) + ".stage2.log")):
        return 1
    for required in [args.output / "mcmctree.out", args.output / "mcmc.txt"]:
        if not required.is_file() or required.stat().st_size == 0:
            raise RuntimeError(f"stage2 missing required output: {required}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
