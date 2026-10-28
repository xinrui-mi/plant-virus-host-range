#!/usr/bin/env bash
set -euo pipefail

FAM="$1"
INPUT="$2"
BASE="$3"
CDHIT="$4"
MAFFT="$5"
TRIMAL="$6"
PYTHON="$7"
THREADS="${8:-12}"

OUTROOT="$BASE/pipeline_results/03K_family_specific_markers/06_AA95_alignment_QC"
D="$OUTROOT/$FAM"
LOG="$D/logs"
mkdir -p "$D" "$LOG"

PREP="$BASE/pipeline/03K_C_prepare_marker_fasta.py"
PARSE="$BASE/pipeline/03K_C_parse_cdhit.py"
FILTER="$BASE/pipeline/03E_filter_alignment_sequences.py"

[[ -s "$INPUT" ]] || { echo "Missing marker input: $INPUT" >&2; exit 10; }
[[ -s "$FILTER" ]] || { echo "Missing frozen alignment filter: $FILTER" >&2; exit 11; }

"$PYTHON" "$PREP" "$INPUT" "$D/${FAM}.marker.clean.faa" \
  > "$LOG/prepare.log" 2>&1

INPUT_N=$(grep -c '^>' "$D/${FAM}.marker.clean.faa" || true)
[[ "$INPUT_N" -ge 10 ]] || { echo "Too few marker sequences: $INPUT_N" >&2; exit 12; }

"$CDHIT" \
  -i "$D/${FAM}.marker.clean.faa" \
  -o "$D/${FAM}.AA95.faa" \
  -c 0.95 \
  -n 5 \
  -d 0 \
  -M 0 \
  -T "$THREADS" \
  > "$LOG/cdhit95.log" 2>&1

REP_N=$(grep -c '^>' "$D/${FAM}.AA95.faa" || true)
[[ "$REP_N" -ge 10 ]] || { echo "Too few AA95 representatives: $REP_N" >&2; exit 13; }

"$PYTHON" "$PARSE" \
  "$D/${FAM}.AA95.faa.clstr" \
  "$D/${FAM}.AA95_membership.tsv" \
  > "$LOG/parse_cdhit.log" 2>&1

"$MAFFT" --auto --thread "$THREADS" \
  "$D/${FAM}.AA95.faa" \
  > "$D/${FAM}.AA95.mafft.faa" \
  2> "$LOG/mafft.log"

ALIGN_N=$(grep -c '^>' "$D/${FAM}.AA95.mafft.faa" || true)
[[ "$ALIGN_N" -ge 10 ]] || { echo "Too few aligned sequences: $ALIGN_N" >&2; exit 14; }

"$TRIMAL" \
  -in "$D/${FAM}.AA95.mafft.faa" \
  -out "$D/${FAM}.AA95.mafft.trim.faa" \
  -automated1 \
  > "$LOG/trimal.log" 2>&1

TRIM_N=$(grep -c '^>' "$D/${FAM}.AA95.mafft.trim.faa" || true)
[[ "$TRIM_N" -ge 10 ]] || { echo "Too few sequences after trimAl: $TRIM_N" >&2; exit 15; }

"$PYTHON" "$FILTER" \
  "$D/${FAM}.AA95.mafft.trim.faa" \
  "$D/${FAM}.AA95.mafft.trim.filtered.faa" \
  --max-gap 0.50 \
  --min-effective-aa 100 \
  --report "$D/${FAM}.alignment_sequence_QC.tsv" \
  > "$LOG/alignment_filter.log" 2>&1

FILTER_N=$(grep -c '^>' "$D/${FAM}.AA95.mafft.trim.filtered.faa" || true)
[[ "$FILTER_N" -ge 10 ]] || { echo "Too few sequences after alignment QC: $FILTER_N" >&2; exit 16; }

{
  echo -e "metric\tvalue"
  echo -e "virus_family\t$FAM"
  echo -e "marker_input_sequence_n\t$INPUT_N"
  echo -e "AA95_representative_n\t$REP_N"
  echo -e "alignment_sequence_n\t$ALIGN_N"
  echo -e "trimmed_sequence_n\t$TRIM_N"
  echo -e "alignment_QC_tree_tip_n\t$FILTER_N"
  echo -e "cdhit_identity\t0.95"
  echo -e "max_gap_fraction\t0.50"
  echo -e "min_effective_aa\t100"
} > "$D/${FAM}.03K_C_summary.tsv"

echo "03K-C PASS $FAM input=$INPUT_N AA95=$REP_N QCtips=$FILTER_N"
