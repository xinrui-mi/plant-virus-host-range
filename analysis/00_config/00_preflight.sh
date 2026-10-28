#!/usr/bin/env bash
# Verifies that the released repository is complete before running any module.
#   bash analysis/00_config/00_preflight.sh
#
# The reproducibility chain starts from the released analysis-ready datasets and
# the released derived inputs shipped under <repo>/data, so no external database
# and no sequence-alignment / phylogeny executable is required for it. The
# executables listed at the end are only needed by the optional pipeline scripts
# that rebuild trees from raw, non-redistributed upstream sequences.
set -euo pipefail
source "$(dirname "$0")/00_paths.sh"
cd "$REPO_ROOT"
Rscript "$ANALYSIS_ROOT/00_config/00_preflight.R"
python3 "$ANALYSIS_ROOT/00_config/00_python_preflight.py"

echo
echo "Optional executables (phylogeny / selection rebuilding only; NOT required for the released chain):"
for x in cd-hit mafft trimal iqtree3 iqtree hmmsearch hmmfetch hyphy HYPHYMPI mpirun; do
  printf '  %-10s %s\n' "$x" "$(command -v "$x" 2>/dev/null || echo 'not found')"
done
printf '\nPreflight PASS\n'
