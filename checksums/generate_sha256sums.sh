#!/usr/bin/env bash
# Regenerate checksums/SHA256SUMS.txt.
#
# Scope: every file that git actually tracks in this repository, which is the
# same set a cloner receives. Anything git-ignored (analysis/*/tables/,
# analysis/*/work/, rendered figure output, *.log, *.contree, ...) is excluded on
# purpose: those are regenerated run artefacts, and a checksum over a file the
# release does not ship would be noise.
#
# Run this LAST, after every other file has stopped changing.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$REPO_ROOT/checksums/SHA256SUMS.txt"

cd "$REPO_ROOT"

if [ ! -d .git ]; then
  echo "ERROR: $REPO_ROOT is not a git work tree; cannot determine the tracked file set." >&2
  exit 1
fi

# Sorted, relative paths, one per line. LC_ALL=C keeps the order stable across
# machines so the file is diffable.
#
# `git ls-files` can emit an entry that is a directory on disk when a nested
# repository or a gitlink is involved; sha256sum refuses those, so they are
# filtered out before hashing rather than being allowed to abort the run.
git ls-files -z \
  | LC_ALL=C sort -z \
  | while IFS= read -r -d '' f; do
      # Skip the checksum file itself: a file cannot contain its own hash.
      [ "$f" = "checksums/SHA256SUMS.txt" ] && continue
      [ -f "$f" ] && printf '%s\0' "$f"
    done \
  | xargs -0 sha256sum -- \
  > "$OUT.tmp"

# sha256sum prefixes each line with "./" only if the path had it; normalise so
# the file always reads "<hash>  <path>" with a plain relative path.
sed -i 's#  \./#  #' "$OUT.tmp"
mv "$OUT.tmp" "$OUT"

printf 'wrote %s: %s files\n' "$OUT" "$(wc -l < "$OUT")"
