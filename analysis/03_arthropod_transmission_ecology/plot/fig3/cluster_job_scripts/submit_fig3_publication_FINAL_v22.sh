#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# ARCHIVAL ONLY — NOT part of the reproduction path.
#
# This is the original LSF submission script from the computing cluster
# (bsub flags, queue names, conda environment, absolute project paths). It is
# kept so that the exact job configuration behind the published figure can be
# inspected, but it will NOT run on another machine. Use the sibling
# `run_*.sh` (or call the plotting R script directly) to regenerate the figure.
# ---------------------------------------------------------------------------
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

bsub \
  -q cpu_bgi \
  -n 1 \
  -R "span[hosts=1]" \
  -J Fig3_pub_v22 \
  -o Fig3_publication_FINAL_v22.%J.out \
  -e Fig3_publication_FINAL_v22.%J.err \
  bash run_fig3_publication_FINAL_v22.sh
