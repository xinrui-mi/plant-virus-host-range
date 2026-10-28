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

PKG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_NAME="Fig2_pub_v14"
LOG_DIR="${PKG_DIR}/bsub_logs"
mkdir -p "${LOG_DIR}"

bsub \
  -q cpu_bgi \
  -n 4 \
  -R "span[hosts=1]" \
  -J "${JOB_NAME}" \
  -o "${LOG_DIR}/${JOB_NAME}.%J.out" \
  -e "${LOG_DIR}/${JOB_NAME}.%J.err" \
  bash "${PKG_DIR}/run_fig2_publication_FINAL_v14.sh"
