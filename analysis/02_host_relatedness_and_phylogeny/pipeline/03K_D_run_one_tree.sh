#!/usr/bin/env bash
set -euo pipefail

FAM="$1"
ALN="$2"
OUTROOT="$3"
IQTREE="$4"
THREADS="${5:-12}"

D="$OUTROOT/$FAM"
mkdir -p "$D/logs"
PRE="$D/$FAM"

[[ -s "$ALN" ]] || { echo "Missing alignment: $ALN" >&2; exit 10; }

N=$(grep -c '^>' "$ALN" || true)
[[ "$N" -ge 10 ]] || { echo "Too few tips: $N" >&2; exit 11; }

"$IQTREE" \
  -s "$ALN" \
  -m MFP \
  -bb 1000 \
  -alrt 1000 \
  -nt "$THREADS" \
  -pre "$PRE" \
  -redo \
  > "$D/logs/iqtree.stdout.log" \
  2> "$D/logs/iqtree.stderr.log"

[[ -s "$PRE.treefile" ]] || { echo "No treefile generated" >&2; exit 12; }
[[ -s "$PRE.iqtree" ]] || { echo "No iqtree report generated" >&2; exit 13; }

echo "03K-D PASS $FAM tips=$N tree=$PRE.treefile"
