#!/usr/bin/env bash
# Regenerate the patch files capturing duckvep's divergence from the pristine
# vendored fastVEP (Huang-lab/fastVEP @ 785922e). One .patch per crate; applying
# them to the upstream src/ reproduces vendor/<crate>/src/. This makes our changes
# reviewable and reproducible (not just prose in ../../PATCHES.md).
#
# Usage: vendor/patches/regen-patches.sh [path-to-fastVEP-checkout-at-785922e]
#   default upstream: a sibling DuckfastVEP checkout at HEAD 785922e
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
VENDOR="$(cd "$HERE/.." && pwd)"
UP="${1:-/root/DuckfastVEP}/crates"
BASE="785922e"

for c in fastvep-core fastvep-genome fastvep-cache fastvep-consequence fastvep-hgvs; do
  [ -d "$UP/$c/src" ] || { echo "skip $c (no upstream src at $UP/$c)"; continue; }
  # a/ = pristine upstream, b/ = our vendored copy — rewrite the absolute paths to
  # a/<crate>/src and b/<crate>/src so `git apply -p1` works from the repo root.
  diff -ruN "$UP/$c/src" "$VENDOR/$c/src" \
    | sed -e "s#${UP}/${c}/src#a/${c}/src#g" -e "s#${VENDOR}/${c}/src#b/${c}/src#g" \
    > "$HERE/$c.patch" || true   # diff exits 1 when files differ
  echo "  $(basename "$HERE/$c.patch")  ($(grep -c '^+' "$HERE/$c.patch" 2>/dev/null || echo 0) added lines)"
done
echo "regenerated patches vs fastVEP@$BASE"
