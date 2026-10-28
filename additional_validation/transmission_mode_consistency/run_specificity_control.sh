#!/usr/bin/env bash
# Specificity control for the transmission-lineage aggregation (NOT_IN_PAPER).
#
# Usage (from anywhere):
#   bash run_specificity_control.sh [analysis_root] [vipdb_input]
#
# Defaults:
#   analysis_root = the `analysis/` directory of this repository
#   vipdb_input   = <repo>/data/analysis_data/VIPDB_relation_normalized.tsv.gz
#
# The original third-party workbook also works as the second argument, e.g.:
#   bash run_specificity_control.sh ../../analysis /path/to/VIPDB_relation_raw.xlsx
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

ANALYSIS_ROOT="${1:-$REPO_ROOT/analysis}"
VIPDB_INPUT="${2:-$REPO_ROOT/data/analysis_data/VIPDB_relation_normalized.tsv.gz}"

R_BIN="${R_BIN:-$(command -v Rscript || true)}"
if [[ -z "$R_BIN" ]]; then
  echo "Rscript not found on PATH. Set R_BIN=/path/to/Rscript." >&2
  exit 1
fi

if [[ ! -f "$VIPDB_INPUT" ]]; then
  echo "Input not found: $VIPDB_INPUT" >&2
  echo "Pass the VIP-DB workbook or the released normalized TSV as the second argument." >&2
  exit 1
fi

exec "$R_BIN" "$SCRIPT_DIR/04D_transmission_specificity_control.R" "$ANALYSIS_ROOT" "$VIPDB_INPUT"
