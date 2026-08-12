#!/usr/bin/env python3
"""Shared launcher for v1 legacy script names."""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


def run(target: list[str], argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    launcher = root / "biohub-rs" / "run-biohub.sh"
    if launcher.exists():
        command = [str(launcher), *target]
    else:
        binary = shutil.which("biohub")
        if binary is None:
            print("biohub not found; build biohub-rs or install biohub on PATH", file=sys.stderr)
            return 127
        command = [binary, *target]
    return subprocess.call([*command, *(sys.argv[1:] if argv is None else argv)])
