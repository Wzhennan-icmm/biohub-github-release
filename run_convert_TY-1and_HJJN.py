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
    cmd = [str(launcher), "gff", "convert-ty1-hjjn", *sys.argv[1:]]
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())
