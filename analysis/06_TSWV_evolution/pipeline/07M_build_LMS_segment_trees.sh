#!/usr/bin/env bash
set -euo pipefail
source "${ANALYSIS_ROOT}/00_config/00_paths.sh"
BASE="${ANALYSIS_ROOT}/06_TSWV_evolution"
D="${BASE}/pipeline_results/06_segment_concordance"
IN="${D}/inputs"; LOG="${BASE}/logs"
mkdir -p "$D" "$IN" "$LOG"
python3 "${BASE}/pipeline/07M_prepare_LMS_segment_inputs.py"

printf 'segment\tmatched_isolate_n\talignment_nt_length\ttreefile\n' > "$D/07M_LMS_segment_tree_manifest.tsv"
for seg in L M S; do
  SD="$D/$seg"; rm -rf "$SD"; mkdir -p "$SD"
  RAW="$IN/${seg}.matched_oriented.fna"
  [[ -s "$RAW" ]] || { echo "Missing $RAW" >&2; exit 20; }
  N=$(grep -c '^>' "$RAW" || true); [[ "$N" -ge 50 ]] || { echo "$seg has <50 matched isolates" >&2; exit 21; }
  mafft --auto --thread "$NCORES" "$RAW" > "$SD/${seg}.mafft.fna" 2> "$LOG/07M_${seg}.mafft.log"
  trimal -in "$SD/${seg}.mafft.fna" -out "$SD/${seg}.trim.fna" -automated1 > "$LOG/07M_${seg}.trimal.log" 2>&1
  N2=$(grep -c '^>' "$SD/${seg}.trim.fna" || true); [[ "$N2" -eq "$N" ]] || { echo "$seg trimAl changed sequence count $N -> $N2" >&2; exit 22; }
  ALEN=$(awk 'BEGIN{s=""}/^>/{if(length(s)>0){print length(s);exit}next}{s=s$0}END{if(length(s)>0)print length(s)}' "$SD/${seg}.trim.fna" | head -1)
  iqtree3 -s "$SD/${seg}.trim.fna" -m MFP -B 1000 --alrt 1000 -bnni -T "$NCORES" --seed 20260813 \
    --prefix "$SD/${seg}" -redo > "$LOG/07M_${seg}.iqtree.stdout.log" 2> "$LOG/07M_${seg}.iqtree.stderr.log"
  [[ -s "$SD/${seg}.treefile" ]] || { echo "$seg treefile missing" >&2; exit 23; }
  printf '%s\t%s\t%s\t%s\n' "$seg" "$N" "$ALEN" "$SD/${seg}.treefile" >> "$D/07M_LMS_segment_tree_manifest.tsv"
done
Rscript "$BASE/analysis/07N_segment_concordance_reassortment.R"
echo "TSWV LMS segment trees + concordance PASS"
