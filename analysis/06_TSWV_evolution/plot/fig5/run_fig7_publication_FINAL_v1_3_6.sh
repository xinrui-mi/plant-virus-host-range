#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# On the cluster this package lived at
# <project>/three_source_manuscript_analysis_v1_1/Fig7_TSWV_case_study/Fig7_publication_FINAL_v1_3_6/.
# In the release those same files are here, so the default is this directory.
PKG_DIR="${PKG_DIR:-$SCRIPT_DIR}"
export FIG_FONT="${FIG_FONT:-Open Sans}"

# Optional: point this at an extra R library if your packages live outside the
# default .libPaths(). Left untouched when unset — reading it as ${VAR} under
# `set -u` would abort the run instead of skipping the setting.
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER
fi

if [[ -z "${RSCRIPT:-}" ]]; then
  RSCRIPT="$(command -v Rscript || true)"
fi
if [[ -z "${RSCRIPT}" ]]; then
  echo "ERROR: Rscript not found. Install R or add it to PATH first." >&2
  exit 1
fi

echo "===== Fig. 5 (TSWV case study) — publication plot ====="
echo "PKG_DIR:        ${PKG_DIR}"
echo "R_LIBS_USER:    ${R_LIBS_USER:-<unset>}"
echo "Rscript:        ${RSCRIPT}"
echo "Font requested: ${FIG_FONT}"
echo
echo "Plotting only. The TSWV analyses and treefiles are read-only frozen inputs."

cd "${PKG_DIR}"

"${RSCRIPT}" 00_check_packages.R
"${RSCRIPT}" 07Q_Fig7_publication_plot_FINAL_v1_3_6.R

echo
echo "Fig. 5 plotting complete."
