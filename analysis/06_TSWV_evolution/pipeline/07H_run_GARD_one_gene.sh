#!/usr/bin/env bash
set -euo pipefail
source "${ANALYSIS_ROOT}/00_config/00_paths.sh"
GENE="${1:?Usage: 07H_run_GARD_one_gene.sh GENE}"
BASE="${ANALYSIS_ROOT}/06_TSWV_evolution"
IN="${BASE}/pipeline_results/03_alignment_and_phylogeny/${GENE}/${GENE}.codon.trim.fna"
D="${BASE}/pipeline_results/04_GARD/${GENE}"
mkdir -p "$D" "${BASE}/logs"
SUB="$D/${GENE}.GARD_input.fna"; AUD="$D/${GENE}.GARD_subset_audit.tsv"
OUT_JSON="$D/${GENE}.GARD.json"; OUT_LF="$D/${GENE}.best-gard.nex"

# Never silently destroy a completed/checkpointed GARD run. Use FORCE_GARD_RERUN=1 only
# when a deliberate clean rerun is intended.
if [[ ( -s "$OUT_JSON" || -s "$OUT_LF" ) && "${FORCE_GARD_RERUN:-0}" != "1" ]]; then
  echo "Existing GARD output detected for $GENE; validating instead of overwriting."
  python3 "$BASE/analysis/07O_finalize_GARD_status.py" --gene "$GENE" --no-write
  exit $?
fi
if [[ "${FORCE_GARD_RERUN:-0}" == "1" ]]; then
  rm -f "$OUT_JSON" "$OUT_LF"
fi

python3 "$BASE/pipeline/07G_make_GARD_diversity_subset.py" "$IN" "$SUB" "$AUD" "${TSWV_GARD_MAX_N:-180}"
[[ -s "$SUB" ]] || { echo "GARD subset FASTA is empty: $SUB" >&2; exit 30; }

# GARD minPartitionSize = 2*N-3. A k-breakpoint model needs k+1 valid partitions.
# Cap the requested breakpoint search to what the alignment can mathematically support;
# this prevents the upstream GARD GA initializer from endlessly proposing impossible models.
read -r NSEQ NSITE MINPART MAXFEAS <<<"$(python3 - "$SUB" <<'PY'
import sys
p=sys.argv[1]; n=0; seqs=[]; s=[]
for line in open(p):
    line=line.strip()
    if line.startswith('>'):
        if s: seqs.append(''.join(s)); s=[]
        n+=1
    elif line: s.append(line)
if s: seqs.append(''.join(s))
L=len(seqs[0]); assert all(len(x)==L for x in seqs)
mp=2*n-3
mf=max(1,L//mp-1)
print(n,L,mp,mf)
PY
)"
REQUESTED="${TSWV_GARD_MAX_BREAKPOINTS:-20}"
if (( REQUESTED < MAXFEAS )); then MAXBP="$REQUESTED"; else MAXBP="$MAXFEAS"; fi
printf 'GARD feasibility %s: n=%s nt=%s minPartition=%s requestedMaxBP=%s effectiveMaxBP=%s\n' "$GENE" "$NSEQ" "$NSITE" "$MINPART" "$REQUESTED" "$MAXBP"

if command -v HYPHYMPI >/dev/null 2>&1 && command -v mpirun >/dev/null 2>&1; then
  RUN=("$(command -v mpirun)" -np "${NCORES:-12}" "$(command -v HYPHYMPI)" gard)
elif [[ "${ALLOW_SERIAL_GARD:-0}" == "1" ]]; then
  RUN=(hyphy gard)
else
  echo "HYPHYMPI and/or mpirun not found" >&2; exit 31
fi

# Nucleotide GARD is used as a recombination screen. Any inferred nucleotide breakpoint
# is subsequently snapped by 07P to the nearest codon boundary (<=1 nt shift) before
# codon-based BUSTED/MEME/FUBAR, while retaining the corresponding GARD partition tree.
"${RUN[@]}" \
  --alignment "$SUB" \
  --type nucleotide \
  --rv GDD \
  --rate-classes 3 \
  --mode Faster \
  --max-breakpoints "$MAXBP" \
  --output "$OUT_JSON" \
  --output-lf "$OUT_LF" \
  'ENV=TOLERATE_NUMERICAL_ERRORS=1;' \
  > "$BASE/logs/${GENE}.GARD.stdout.md" \
  2> "$BASE/logs/${GENE}.GARD.stderr.log"

[[ -s "$OUT_JSON" ]] || { echo "GARD JSON missing/empty: $OUT_JSON" >&2; exit 32; }
[[ -s "$OUT_LF" ]] || { echo "GARD snapshot missing/empty: $OUT_LF" >&2; exit 33; }
python3 "$BASE/analysis/07O_finalize_GARD_status.py" --gene "$GENE" --no-write

echo "GARD scientific gate PASS $GENE"
