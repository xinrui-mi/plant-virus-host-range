#!/usr/bin/env bash
set -euo pipefail

fam="$1"
pfam="$2"
hmmscan_bin="$3"
base="$4"
ncores="${5:-12}"

ORF="$base/pipeline_results/03K_family_specific_markers/01_orf_recovery"
OUT="$base/pipeline_results/03K_family_specific_markers/02_pfam"
LOG="$base/pipeline_results/03K_family_specific_markers/logs"

FAA="$ORF/${fam}.candidate_orfs.faa"
DOM="$OUT/${fam}.Pfam.domtblout"

[[ -s "$FAA" ]] || {
  echo "ERROR: missing ORF FASTA: $FAA" >&2
  exit 11
}

"$hmmscan_bin" \
  --cpu "$ncores" \
  --noali \
  --domtblout "$DOM" \
  "$pfam" \
  "$FAA" \
  > "$LOG/${fam}.Pfam.hmmscan.stdout.log" \
  2> "$LOG/${fam}.Pfam.hmmscan.stderr.log"

[[ -s "$DOM" ]] || {
  echo "ERROR: hmmscan produced empty/missing domtblout: $DOM" >&2
  exit 12
}

echo "PASS $fam"
echo "PFAM_HIT_ROWS=$(grep -vc '^#' "$DOM" || true)"
