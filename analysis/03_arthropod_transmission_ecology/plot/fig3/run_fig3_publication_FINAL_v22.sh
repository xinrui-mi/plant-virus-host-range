#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# On the cluster this directory sat inside a hand-built tree
# (<project>/three_source_manuscript_analysis_v1_1/...). In the release the same
# module lives at <repo>/analysis/<module>/, so the defaults point there. Both
# are overridable, which is the only thing that ever differed between them.
export FINAL_THREE_SOURCE_ROOT="${FINAL_THREE_SOURCE_ROOT:-$(cd "$SCRIPT_DIR/../../../.." && pwd)}"
export FIG3_ANALYSIS_ROOT="${FIG3_ANALYSIS_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
export FINAL_PIPELINE_RESULTS="${FINAL_PIPELINE_RESULTS:-$FINAL_THREE_SOURCE_ROOT/data/analysis_data}"
export FIG_FONT="${FIG_FONT:-Open Sans}"

# Optional: point this at an extra R library if your packages live outside the
# default .libPaths(). Left untouched when unset — testing for it with ${VAR}
# under `set -u` would abort the run instead of skipping the setting.
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER
fi

RSCRIPT="${RSCRIPT:-$(command -v Rscript)}"

echo "FINAL_THREE_SOURCE_ROOT: $FINAL_THREE_SOURCE_ROOT"
echo "FINAL_PIPELINE_RESULTS: $FINAL_PIPELINE_RESULTS"
echo "FIG3_ANALYSIS_ROOT: $FIG3_ANALYSIS_ROOT"
echo "R_LIBS_USER: ${R_LIBS_USER:-<unset>}"
echo "Rscript: $RSCRIPT"
echo "Font requested: $FIG_FONT"

"$RSCRIPT" "$SCRIPT_DIR/03_Fig3_publication_plot_FINAL_v22.R"
