#!/usr/bin/env bash
set -euo pipefail
source "${ANALYSIS_ROOT}/00_config/00_paths.sh"
GENE="${1:?Usage: 07I_run_selection_one_gene.sh GENE}"
BASE="${ANALYSIS_ROOT}/06_TSWV_evolution"
STATUS="${BASE}/tables/16_TSWV_GARD_final_status.tsv"
INPUT="${BASE}/pipeline_results/05_selection/00_inputs/${GENE}/${GENE}.selection_codon_safe.nex"
OUT="${BASE}/pipeline_results/05_selection/${GENE}"
mkdir -p "$OUT" "${BASE}/logs"

[[ -s "$STATUS" ]] || { echo "Missing GARD status table: $STATUS" >&2; exit 40; }
[[ -s "$INPUT" ]] || { echo "Missing codon-safe selection input: $INPUT" >&2; exit 41; }

python3 - "$STATUS" "$GENE" <<'PY'
import csv,sys
p,g=sys.argv[1:]
valid={'FINAL_GARD','PASS_NO_BREAKPOINT_MAX_FEASIBLE_BP1','PASS_SINGLE_BREAKPOINT_MAX_FEASIBLE_BP1'}
with open(p,newline='') as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
r=next((x for x in rows if x['gene']==g),None)
if not r: raise SystemExit(f'No GARD status row for {g}')
if r['status'] not in valid: raise SystemExit(f'Unsafe GARD status for {g}: {r["status"]}')
print('Selection gate PASS',g,r['status'],'breakpoints=',r.get('final_breakpoint_n'))
PY

rm -f "$OUT/${GENE}.BUSTED.json" "$OUT/${GENE}.MEME.json" "$OUT/${GENE}.FUBAR.json"

hyphy busted --alignment "$INPUT" --branches All --output "$OUT/${GENE}.BUSTED.json" \
  > "$BASE/logs/${GENE}.BUSTED.stdout.md" 2> "$BASE/logs/${GENE}.BUSTED.stderr.log"

hyphy meme --alignment "$INPUT" --branches All --output "$OUT/${GENE}.MEME.json" \
  > "$BASE/logs/${GENE}.MEME.stdout.md" 2> "$BASE/logs/${GENE}.MEME.stderr.log"

hyphy fubar --alignment "$INPUT" --output "$OUT/${GENE}.FUBAR.json" \
  > "$BASE/logs/${GENE}.FUBAR.stdout.md" 2> "$BASE/logs/${GENE}.FUBAR.stderr.log"

python3 - "$OUT/${GENE}.BUSTED.json" "$OUT/${GENE}.MEME.json" "$OUT/${GENE}.FUBAR.json" <<'PY'
import json,sys,os
for p in sys.argv[1:]:
    if not os.path.exists(p) or os.path.getsize(p)==0:
        raise SystemExit('Missing/empty selection JSON: '+p)
    with open(p) as f: json.load(f)
print('Selection JSON validation PASS')
PY

echo "Selection PASS $GENE"
