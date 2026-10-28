#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Fig. 4 reads two frozen analysis modules: the network itself (04) and the
# VIP-DB vector adjustment panel, which lives in the candidate-contextualisation
# module (05). On the cluster those were Fig5_network_prioritization and
# Fig6_candidate_contextualization; here they are analysis/04_ and analysis/05_.
# Earlier revisions of this wrapper pointed CONTEXT_ROOT at the network module,
# so the adjustment table was looked for one directory away from where it is.
NETWORK_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONTEXT_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/../../../05_vector_contextualization" 2>/dev/null && pwd || echo "")"
ASSOC_MASTER_DEFAULT="$(cd "$SCRIPT_DIR/../../../.." && pwd)/data/analysis_input/association_analysis_input.tsv.gz"

export NETWORK_ROOT="${NETWORK_ROOT:-$NETWORK_ROOT_DEFAULT}"
export CONTEXT_ROOT="${CONTEXT_ROOT:-$CONTEXT_ROOT_DEFAULT}"
export ASSOC_MASTER="${ASSOC_MASTER:-$ASSOC_MASTER_DEFAULT}"
# Optional: point this at an extra R library if your packages live outside the
# default .libPaths(). Left untouched when unset — reading it as ${VAR} under
# `set -u` would abort the run instead of skipping the setting.
if [ -n "${R_LIBS_USER:-}" ]; then export R_LIBS_USER; fi
export FIG_FONT="${FIG_FONT:-Open Sans}"

if [[ -z "${RSCRIPT:-}" ]]; then
  RSCRIPT="$(command -v Rscript || true)"
fi
if [[ -z "${RSCRIPT}" ]]; then echo "ERROR: Rscript not found. Install R or add it to PATH first." >&2; exit 1; fi

echo "NETWORK_ROOT: ${NETWORK_ROOT}"
echo "CONTEXT_ROOT: ${CONTEXT_ROOT}"
echo "ASSOC_MASTER: ${ASSOC_MASTER}"
echo "R_LIBS_USER:  ${R_LIBS_USER:-<unset>}"
echo "Rscript:      ${RSCRIPT}"
"${RSCRIPT}" "${SCRIPT_DIR}/04_Fig4_publication_FINAL_v1_1.R"
