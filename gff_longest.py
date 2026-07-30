#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    launcher = Path(__file__).resolve().parent / "biohub-rs" / "run-biohub.sh"
    if not launcher.exists():
        print("biohub-rs launcher not found: {}".format(launcher), file=sys.stderr)
        return 127

    if len(sys.argv) == 1:
        print("Usage: gff_longest.py <input_gff> [output_gff]")
        return 1
    if len(sys.argv) >= 2 and len(sys.argv) <= 3:
        argv = [str(launcher), "gff", "filter-gemoma"]
        argv += ["--gff", sys.argv[1]]
        if len(sys.argv) == 3:
            argv += ["--output", sys.argv[2]]
        return subprocess.call(argv)
    print("Usage: gff_longest.py <input_gff> [output_gff]")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
