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

PROJECT_ROOT="/share/home/bgi_mixr/project/biye/0_newstart/0_rawseq/final_three_source"
ANALYSIS_ROOT="${PROJECT_ROOT}/three_source_manuscript_analysis_v1_1"
FIG7_ROOT="${ANALYSIS_ROOT}/Fig7_TSWV_case_study"
PKG_DIR="${FIG7_ROOT}/Fig7_publication_FINAL_v1_3_6"
LOG_DIR="${FIG7_ROOT}/publication_figures_FINAL_v1_3_6/logs"

mkdir -p "${LOG_DIR}"

bsub \
  -q cpu_bgi \
  -n 1 \
  -R "span[hosts=1]" \
  -J Fig7Plot_v1_3_6 \
  -o "${LOG_DIR}/Fig7Plot_v1_3_6.%J.out" \
  -e "${LOG_DIR}/Fig7Plot_v1_3_6.%J.err" \
  -cwd "${PKG_DIR}" \
  bash "${PKG_DIR}/run_fig7_publication_FINAL_v1_3_6.sh"
