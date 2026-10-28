#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
BASE="$ROOT/three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint"
KROOT="$BASE/pipeline_results/03K_family_specific_markers"
OUT="$KROOT/07_Geminiviridae_Rep_conserved"
D="$OUT/AA95_alignment_QC"
mkdir -p "$D/logs"

python3 "$BASE/pipeline/03K_C3_extract_Geminiviridae_Rep_conserved.py" \
  --project-root "$ROOT" | tee "$OUT/extraction.stdout.log"

INPUT="$OUT/Geminiviridae.Rep_conserved.raw.faa"
[[ -s "$INPUT" ]] || { echo "Missing extracted conserved Rep marker" >&2; exit 2; }

# Freeze data-rich criterion before tree building.
SPECIES_N=$(awk -F'\t' '$1=="species_level_label_n"{print $2}' "$OUT/Geminiviridae.Rep_conserved.summary.tsv")
if [[ "$SPECIES_N" -lt 100 ]]; then
  echo "FAIL: conserved Rep marker covers fewer than 100 species-level labels: $SPECIES_N" >&2
  exit 3
fi

python3 "$BASE/pipeline/03K_C_prepare_marker_fasta.py" \
  "$INPUT" "$D/Geminiviridae.Rep_conserved.clean.faa" \
  > "$D/logs/prepare.log" 2>&1

cd-hit \
  -i "$D/Geminiviridae.Rep_conserved.clean.faa" \
  -o "$D/Geminiviridae.Rep_conserved.AA95.faa" \
  -c 0.95 -n 5 -d 0 -M 0 -T 12 \
  > "$D/logs/cdhit95.log" 2>&1

python3 "$BASE/pipeline/03K_C_parse_cdhit.py" \
  "$D/Geminiviridae.Rep_conserved.AA95.faa.clstr" \
  "$D/Geminiviridae.Rep_conserved.AA95_membership.tsv" \
  > "$D/logs/parse_cdhit.log" 2>&1

mafft --auto --thread 12 \
  "$D/Geminiviridae.Rep_conserved.AA95.faa" \
  > "$D/Geminiviridae.Rep_conserved.AA95.mafft.faa" \
  2> "$D/logs/mafft.log"

trimal \
  -in "$D/Geminiviridae.Rep_conserved.AA95.mafft.faa" \
  -out "$D/Geminiviridae.Rep_conserved.AA95.mafft.trim.faa" \
  -automated1 \
  > "$D/logs/trimal.log" 2>&1

python3 "$BASE/pipeline/03E_filter_alignment_sequences.py" \
  "$D/Geminiviridae.Rep_conserved.AA95.mafft.trim.faa" \
  "$D/Geminiviridae.Rep_conserved.AA95.mafft.trim.filtered.faa" \
  --max-gap 0.50 \
  --min-effective-aa 100 \
  --report "$D/Geminiviridae.Rep_conserved.alignment_sequence_QC.tsv" \
  > "$D/logs/alignment_filter.log" 2>&1

IN_N=$(grep -c '^>' "$INPUT" || true)
AA95_N=$(grep -c '^>' "$D/Geminiviridae.Rep_conserved.AA95.faa" || true)
QC_N=$(grep -c '^>' "$D/Geminiviridae.Rep_conserved.AA95.mafft.trim.filtered.faa" || true)

{
  echo -e "metric\tvalue"
  echo -e "marker_input_sequence_n\t$IN_N"
  echo -e "species_level_label_n\t$SPECIES_N"
  echo -e "AA95_representative_n\t$AA95_N"
  echo -e "alignment_QC_tree_tip_n\t$QC_N"
  echo -e "cdhit_identity\t0.95"
  echo -e "max_gap_fraction\t0.50"
  echo -e "min_effective_aa\t100"
} > "$OUT/Geminiviridae.Rep_conserved.03K_C3_summary.tsv"

echo "03K-C3 PASS input=$IN_N species=$SPECIES_N AA95=$AA95_N QCtips=$QC_N"
