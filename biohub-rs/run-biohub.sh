#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPECTED_VERSION="$(awk -F '"' '/^version =/ {print $2; exit}' "$ROOT/Cargo.toml")"

is_current_binary() {
  local candidate="$1"
  [ -x "$candidate" ] || return 1
  [ "$candidate" -nt "$ROOT/Cargo.toml" ] || return 1
  [ "$candidate" -nt "$ROOT/Cargo.lock" ] || return 1
  if find "$ROOT/src" -type f -newer "$candidate" -print -quit | grep -q .; then
    return 1
  fi
  [ "$("$candidate" --version)" = "$EXPECTED_VERSION" ]
}

for candidate in \
  "$ROOT/target/release/biohub" \
  "$ROOT/target/debug/biohub" \
  "$ROOT/target/biohub" \
  "$ROOT/target/release/biohub-rs" \
  "$ROOT/target/debug/biohub-rs"
do
  if is_current_binary "$candidate"; then
    exec "$candidate" "$@"
  fi
done

if command -v cargo >/dev/null 2>&1; then
  echo "Building BioHub $EXPECTED_VERSION from current source..." >&2
  cargo build --release --locked --manifest-path "$ROOT/Cargo.toml"
  exec "$ROOT/target/release/biohub" "$@"
fi

echo "Current BioHub $EXPECTED_VERSION binary not found and cargo is unavailable." >&2
exit 127
