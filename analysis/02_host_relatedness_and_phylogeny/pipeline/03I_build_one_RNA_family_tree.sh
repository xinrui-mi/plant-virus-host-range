#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <virus_family>" >&2
  exit 2
fi

FAM="$1"

FAMILIES=(
  Potyviridae
  Betaflexiviridae
  Secoviridae
  Tombusviridae
  Rhabdoviridae
  Solemoviridae
  Partitiviridae
  Closteroviridae
  Mitoviridae
  Alphaflexiviridae
)

valid=0
for x in "${FAMILIES[@]}"; do
  [[ "$FAM" == "$x" ]] && valid=1
done
[[ "$valid" -eq 1 ]] || { echo "Unsupported family: $FAM" >&2; exit 2; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT="$(cd "${BASE}/../.." && pwd)"

WORK="${BASE}/work/03I_RNA_standardized_RdRP"
OUTROOT="${BASE}/pipeline_results/03I_RNA_standardized_RdRP"
D="${OUTROOT}/${FAM}"
LOG="${BASE}/logs/03I_RNA_standardized_RdRP"
TABLE="${BASE}/tables"

mkdir -p "$D" "$LOG" "$TABLE"

RAW="${WORK}/${FAM}.RdRP_core.raw.faa"
[[ -s "$RAW" ]] || {
  echo "Missing/empty input: $RAW" >&2
  echo "Run 03I_prepare_RNA_family_RdRP_core.py first." >&2
  exit 3
}

NCORES="${FIG3_03I_NCORES:-12}"
MEM_MB="${FIG3_03I_CDHIT_MEMORY_MB:-16000}"
SEED="${FIG3_03I_SEED:-20260922}"

command -v cd-hit >/dev/null 2>&1 || { echo "cd-hit not found" >&2; exit 4; }
command -v mafft >/dev/null 2>&1 || { echo "mafft not found" >&2; exit 4; }
command -v trimal >/dev/null 2>&1 || { echo "trimal not found" >&2; exit 4; }

if command -v iqtree3 >/dev/null 2>&1; then
  IQTREE="$(command -v iqtree3)"
elif command -v iqtree2 >/dev/null 2>&1; then
  IQTREE="$(command -v iqtree2)"
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE="$(command -v iqtree)"
else
  echo "IQ-TREE executable not found (iqtree3/iqtree2/iqtree)" >&2
  exit 4
fi

FILTER="${BASE}/pipeline/03E_filter_alignment_sequences.py"
PARSER="${BASE}/pipeline/03I_parse_cdhit.py"
[[ -s "$FILTER" ]] || { echo "Missing frozen alignment filter: $FILTER" >&2; exit 4; }
[[ -s "$PARSER" ]] || { echo "Missing 03I CD-HIT parser: $PARSER" >&2; exit 4; }

echo "===== 03I ${FAM} ====="
echo "Input: $RAW"
echo "Cores: $NCORES"
echo "IQ-TREE: $IQTREE"

RAW_N=$(grep -c '^>' "$RAW" || true)
[[ "$RAW_N" -ge 2 ]] || { echo "Too few raw marker sequences: $RAW_N" >&2; exit 5; }

# ------------------------------------------------------------------
# 1. 95% AA clustering — same analysis concept as old 03E.
#    No upstream/downstream extension; the input itself is the frozen
#    standardized Stage-8 RdRP conserved core.
# ------------------------------------------------------------------
NR95="${D}/${FAM}.RdRP_core.nr95.faa"

if [[ "${FIG3_03I_REUSE:-1}" == "1" && -s "$NR95" && -s "${NR95}.clstr" ]]; then
  echo "[REUSE] CD-HIT output"
else
  rm -f "$NR95" "${NR95}.clstr"
  cd-hit \
    -i "$RAW" \
    -o "$NR95" \
    -c 0.95 \
    -n 5 \
    -d 0 \
    -g 1 \
    -T "$NCORES" \
    -M "$MEM_MB" \
    > "${LOG}/${FAM}.cdhit95.log" 2>&1
fi

python3 "$PARSER" \
  --family "$FAM" \
  --clstr "${NR95}.clstr" \
  --out "${D}/${FAM}.cluster95_membership.tsv" \
  > "${LOG}/${FAM}.cdhit95.parse.log" 2>&1

NR95_N=$(grep -c '^>' "$NR95" || true)
[[ "$NR95_N" -ge 2 ]] || { echo "Too few NR95 sequences: $NR95_N" >&2; exit 6; }

# ------------------------------------------------------------------
# 2. MAFFT
# ------------------------------------------------------------------
MAFFT_OUT="${D}/${FAM}.RdRP_core.nr95.mafft.faa"
if [[ "${FIG3_03I_REUSE:-1}" == "1" && -s "$MAFFT_OUT" ]]; then
  echo "[REUSE] MAFFT"
else
  mafft --auto --thread "$NCORES" "$NR95" \
    > "$MAFFT_OUT" \
    2> "${LOG}/${FAM}.mafft.log"
fi

# ------------------------------------------------------------------
# 3. trimAl automated1
# ------------------------------------------------------------------
TRIM="${D}/${FAM}.RdRP_core.nr95.trim.faa"
if [[ "${FIG3_03I_REUSE:-1}" == "1" && -s "$TRIM" ]]; then
  echo "[REUSE] trimAl"
else
  trimal \
    -in "$MAFFT_OUT" \
    -out "$TRIM" \
    -automated1 \
    > "${LOG}/${FAM}.trimal.log" 2>&1
fi

# ------------------------------------------------------------------
# 4. Same sequence-level alignment QC as frozen 03E
#    max gap fraction <= 0.50; effective length >= 100 aa
# ------------------------------------------------------------------
FILTERED="${D}/${FAM}.RdRP_core.nr95.trim.filtered.faa"
QCSEQ="${D}/${FAM}.alignment_sequence_QC.tsv"

python3 "$FILTER" \
  "$TRIM" \
  "$FILTERED" \
  --max-gap "${FIG3_TREE_MAX_GAP_FRACTION:-0.50}" \
  --min-effective-aa "${FIG3_TREE_MIN_EFFECTIVE_AA:-100}" \
  --report "$QCSEQ" \
  > "${LOG}/${FAM}.alignment_filter.log" 2>&1

NF=$(grep -c '^>' "$FILTERED" || true)
[[ "$NF" -ge 4 ]] || {
  echo "Too few sequences after alignment QC for a phylogeny: $NF" >&2
  exit 7
}

# ------------------------------------------------------------------
# 5. IQ-TREE + ModelFinder + UFBoot1000 + SH-aLRT1000
# ------------------------------------------------------------------
PREFIX="${D}/${FAM}.RdRP_core"
TREEFILE="${PREFIX}.treefile"

if [[ "${FIG3_03I_REUSE:-1}" == "1" && -s "$TREEFILE" ]]; then
  echo "[REUSE] IQ-TREE treefile"
else
  "$IQTREE" \
    -s "$FILTERED" \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    -T "$NCORES" \
    -seed "$SEED" \
    --prefix "$PREFIX" \
    --redo \
    > "${LOG}/${FAM}.iqtree.stdout.log" \
    2> "${LOG}/${FAM}.iqtree.stderr.log"
fi

[[ -s "$TREEFILE" ]] || { echo "Treefile missing: $TREEFILE" >&2; exit 8; }

# ------------------------------------------------------------------
# 6. Per-family manifest
# ------------------------------------------------------------------
MANIFEST="${D}/${FAM}.03I_manifest.tsv"
{
  printf "virus_family\tmarker\traw_core_n\tcluster95_n\tpost_alignment_QC_tip_n\tidentity_threshold\talignment_filter_max_gap\talignment_filter_min_effective_aa\ttreefile\tstatus\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n" \
    "$FAM" \
    "standardized Stage-8 RdRP conserved core" \
    "$RAW_N" \
    "$NR95_N" \
    "$NF" \
    "0.95" \
    "${FIG3_TREE_MAX_GAP_FRACTION:-0.50}" \
    "${FIG3_TREE_MIN_EFFECTIVE_AA:-100}" \
    "$TREEFILE"
} > "$MANIFEST"

echo "03I PASS ${FAM}: raw=${RAW_N} nr95=${NR95_N} filtered_tips=${NF}"
