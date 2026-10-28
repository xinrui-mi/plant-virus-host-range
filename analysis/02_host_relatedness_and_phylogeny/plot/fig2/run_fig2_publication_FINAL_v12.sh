#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIG3_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/../.." && pwd)"
export FIG3_ROOT="${FIG3_ROOT:-$FIG3_ROOT_DEFAULT}"
# Optional: set R_LIBS_USER to an additional R library if your packages live outside the default.
export FIG_FONT="${FIG_FONT:-Open Sans}"

echo "FIG3_ROOT: ${FIG3_ROOT}"
echo "R_LIBS_USER: ${R_LIBS_USER:-<unset>}"
echo "Rscript: $(command -v Rscript)"
echo "Font requested: ${FIG_FONT}"

Rscript "${SCRIPT_DIR}/03_Fig2_publication_plot_FINAL_v3.R"
