#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for candidate in \
  "$ROOT/target/release/biohub" \
  "$ROOT/target/debug/biohub" \
  "$ROOT/target/biohub" \
  "$ROOT/target/release/biohub-rs" \
  "$ROOT/target/debug/biohub-rs"
do
  if [ -x "$candidate" ]; then
    exec "$candidate" "$@"
  fi
done

if command -v biohub >/dev/null 2>&1; then
  exec biohub "$@"
fi

echo "biohub-rs binary not found. Build it first: (cd \"$ROOT\" && cargo build --release)." >&2
exit 127

