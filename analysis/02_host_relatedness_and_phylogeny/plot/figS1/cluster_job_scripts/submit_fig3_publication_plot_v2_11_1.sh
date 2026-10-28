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
FIG3_ROOT_DEFAULT="/share/home/bgi_mixr/project/biye/0_newstart/0_rawseq/final_three_source/three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint"

export FIG3_ROOT="${FIG3_ROOT:-$FIG3_ROOT_DEFAULT}"

mkdir -p "${FIG3_ROOT}/publication_figures_v2_11_1/logs"

bsub \
  -q cpu_bgi \
  -n 2 \
  -R "span[hosts=1]" \
  -J Fig3_plot_v2_11_1 \
  -o "${FIG3_ROOT}/publication_figures_v2_11_1/logs/Fig3_plot_v2_11_1.%J.out" \
  -e "${FIG3_ROOT}/publication_figures_v2_11_1/logs/Fig3_plot_v2_11_1.%J.err" \
  -cwd "${SCRIPT_DIR}" \
  bash "${SCRIPT_DIR}/run_fig3_publication_plot_v2_11_1.sh"
