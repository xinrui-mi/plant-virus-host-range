#!/usr/bin/env bash
# Standalone Supplementary Fig. S1 wrapper.
# Draws ONLY the submitted S1 (four-panel composite A/B/C/D) with the frozen
# pre-final logic extracted from 03_Fig3_publication_plot_v2_11_1.R.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

FIG3_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/../.." && pwd)"
export FIG3_ROOT="${FIG3_ROOT:-$FIG3_ROOT_DEFAULT}"

echo "FIG3_ROOT: ${FIG3_ROOT}"
echo "R_LIBS_USER: ${R_LIBS_USER:-<unset>}"
echo "Rscript: $(command -v Rscript)"
export FIG_FONT="${FIG_FONT:-Open Sans}"
echo "Font requested: ${FIG_FONT}"

Rscript "${SCRIPT_DIR}/01_FigS1_publication_plot_v2_11_1_standalone.R"
