#!/usr/bin/env python3
from __future__ import annotations

import sys

from biohub_legacy import run


if len(sys.argv) not in (2, 3):
    print("Usage: gff_longest.py <input_gff> [output_gff]", file=sys.stderr)
    raise SystemExit(2)

args = ["--gff", sys.argv[1]]
if len(sys.argv) == 3:
    args.extend(["--output", sys.argv[2]])
raise SystemExit(run(["gff", "filter-gemoma"], args))
