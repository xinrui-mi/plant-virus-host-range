#!/usr/bin/env bash
set -euo pipefail
source "${ANALYSIS_ROOT}/00_config/00_paths.sh"
BASE="${ANALYSIS_ROOT}/02_host_relatedness_and_phylogeny"
WORK="${BASE}/work/marker_inputs"
OUT="${BASE}/pipeline_results"
LOG="${BASE}/logs"
mkdir -p "$WORK" "$OUT" "$LOG"

RDRPSCAN_HMM="${RDRPSCAN_HMM:-${HOME}/source/VirusMarker_HMM/rdrpscan_full.hmm}"
[[ -s "$PFAM_HMM" ]] || { echo "Missing PFAM_HMM: $PFAM_HMM" >&2; exit 2; }
[[ -s "${PFAM_HMM}.h3m" ]] || { echo "Pfam database is not hmmpress-indexed: ${PFAM_HMM}.h3m missing" >&2; exit 2; }
[[ -s "$RDRPSCAN_HMM" ]] || { echo "Missing RDRPSCAN_HMM: $RDRPSCAN_HMM" >&2; exit 2; }
[[ -s "$ALL_THREE_FASTA" ]] || { echo "Missing ALL_THREE_FASTA: $ALL_THREE_FASTA" >&2; exit 2; }

# Remove only obsolete/failed 03E implementation artefacts; preserve 03A-03D outputs.
rm -rf "$WORK/pfam_marker_models"
rm -f "$WORK/PF00799.hmm" "$WORK/PF00078.hmm" "$WORK/PF00075.hmm" \
      "$WORK/Geminiviridae.PF00799.domtbl" "$WORK/Caulimoviridae.PF00078.domtbl" "$WORK/Caulimoviridae.PF00075.domtbl" \
      "$WORK/Geminiviridae.GEM_REP.domtbl" "$WORK/Caulimoviridae.CAUL_RT.domtbl" "$WORK/Caulimoviridae.CAUL_RH.domtbl"

# -------- DNA conserved replication markers --------
# 1) extract all CDS for every eligible sequence;
# 2) annotation-screen to <=2 biologically plausible proteins per sequence;
# 3) scan those candidates against the complete, version-independent Pfam-A database;
# 4) identify domains semantically from current Pfam NAME/DESC and extract conserved marker regions.
python3 "${BASE}/pipeline/03E_extract_dna_cds.py"
python3 "${BASE}/pipeline/03E_select_dna_candidates.py"

[[ -s "$WORK/Geminiviridae.Rep_candidates.faa" ]] || { echo "Empty Geminiviridae Rep candidate FASTA" >&2; exit 5; }
[[ -s "$WORK/Caulimoviridae.Pol_candidates.faa" ]] || { echo "Empty Caulimoviridae Pol candidate FASTA" >&2; exit 5; }

hmmscan --cpu "$NCORES" --noali --domtblout "$WORK/Geminiviridae.Pfam.domtbl" \
  "$PFAM_HMM" "$WORK/Geminiviridae.Rep_candidates.faa" \
  > "$LOG/Geminiviridae.full_Pfam.hmmscan.log" 2>&1
hmmscan --cpu "$NCORES" --noali --domtblout "$WORK/Caulimoviridae.Pfam.domtbl" \
  "$PFAM_HMM" "$WORK/Caulimoviridae.Pol_candidates.faa" \
  > "$LOG/Caulimoviridae.full_Pfam.hmmscan.log" 2>&1
python3 "${BASE}/pipeline/03E_extract_dna_markers_fullpfam.py"

# -------- RNA RdRP markers --------
# Family-specific windows are retained as pre-defined marker windows, with no old taxon cap.
python3 "${BASE}/pipeline/03E_prepare_rna_orfs.py"
for fam in Potyviridae Tombusviridae; do
  hmmsearch --cpu "$NCORES" --noali --domtblout "$WORK/${fam}.RdRP.domtbl" \
    "$RDRPSCAN_HMM" "$WORK/${fam}.all_ORFs.faa" > "$LOG/${fam}.RdRP.hmmsearch.log" 2>&1
done
python3 "${BASE}/pipeline/03E_extract_rna_markers.py"

declare -A RAW MARKER
RAW[Geminiviridae]="$WORK/Geminiviridae.Rep_core.raw.faa"; MARKER[Geminiviridae]="Rep conserved region (full-Pfam verified)"
RAW[Caulimoviridae]="$WORK/Caulimoviridae.RT_RNaseH.raw.faa"; MARKER[Caulimoviridae]="RT-to-RNaseH conserved region (full-Pfam verified)"
RAW[Potyviridae]="$WORK/Potyviridae.RdRP_expanded.raw.faa"; MARKER[Potyviridae]="RdRP HMM hit + upstream 220 aa / downstream 120 aa"
RAW[Tombusviridae]="$WORK/Tombusviridae.RdRP_expanded.raw.faa"; MARKER[Tombusviridae]="RdRP HMM hit + upstream 80 aa / downstream 180 aa"

printf 'virus_family\tmarker_type\traw_marker_n\tnr95_n\tfiltered_tip_n\talignment_aa_length\n' > "$OUT/03E_tree_pipeline_manifest.tsv"
for fam in Geminiviridae Caulimoviridae Potyviridae Tombusviridae; do
  IN="${RAW[$fam]}"
  [[ -s "$IN" ]] || { echo "$fam marker FASTA missing or empty: $IN" >&2; exit 10; }
  N=$(grep -c '^>' "$IN" || true)
  [[ "$N" -ge 4 ]] || { echo "$fam has <4 marker sequences" >&2; exit 10; }

  # A complete rerun is intentional: remove only this family's 03E outputs so no old checkpoint/tree can leak in.
  D="$OUT/$fam"
  rm -rf "$D"
  mkdir -p "$D"
  awk '/^>/{print;next}{gsub(/[JBOZU]/,"X");gsub(/\*/,"");print}' "$IN" > "$D/${fam}.marker.clean.faa"

  cd-hit -i "$D/${fam}.marker.clean.faa" -o "$D/${fam}.nr95.faa" \
    -c 0.95 -n 5 -G 1 -aS 0.80 -aL 0.80 -g 1 -d 0 -T "$NCORES" -M 0 \
    > "$LOG/${fam}.cdhit.log" 2>&1
  python3 "$BASE/pipeline/03E_parse_cdhit.py" "$D/${fam}.nr95.faa.clstr" "$fam" "$D/${fam}.nr95.membership.tsv"
  NR=$(grep -c '^>' "$D/${fam}.nr95.faa" || true)
  [[ "$NR" -ge 4 ]] || { echo "$fam has <4 CD-HIT95 representatives" >&2; exit 11; }

  if [[ "$NR" -le 300 ]]; then
    mafft --localpair --maxiterate 1000 --thread "$NCORES" "$D/${fam}.nr95.faa" > "$D/${fam}.mafft.faa" 2> "$LOG/${fam}.mafft.log"
  else
    mafft --auto --thread "$NCORES" "$D/${fam}.nr95.faa" > "$D/${fam}.mafft.faa" 2> "$LOG/${fam}.mafft.log"
  fi
  trimal -in "$D/${fam}.mafft.faa" -out "$D/${fam}.trim.faa" -automated1 > "$LOG/${fam}.trimal.log" 2>&1
  python3 "$BASE/pipeline/03E_filter_alignment_sequences.py" "$D/${fam}.trim.faa" "$D/${fam}.trim.filtered.faa" \
    --max-gap "${FIG3_TREE_MAX_GAP_FRACTION:-0.50}" --min-effective-aa "${FIG3_TREE_MIN_EFFECTIVE_AA:-100}" \
    --report "$D/${fam}.alignment_sequence_QC.tsv" > "$LOG/${fam}.alignment_filter.log" 2>&1
  NF=$(grep -c '^>' "$D/${fam}.trim.filtered.faa" || true)
  [[ "$NF" -ge 4 ]] || { echo "$fam has <4 tips after alignment QC" >&2; exit 12; }
  ALEN=$(awk 'BEGIN{s=""}/^>/{if(length(s)>0){print length(s);exit}next}{s=s$0}END{if(length(s)>0)print length(s)}' "$D/${fam}.trim.filtered.faa" | head -1)

  iqtree3 -s "$D/${fam}.trim.filtered.faa" -m MFP -B 1000 --alrt 1000 -bnni \
    -T "$NCORES" --seed 20260813 --prefix "$D/${fam}" -redo \
    > "$LOG/${fam}.iqtree.stdout.log" 2> "$LOG/${fam}.iqtree.stderr.log"
  [[ -s "$D/${fam}.treefile" ]] || { echo "$fam IQ-TREE treefile missing" >&2; exit 13; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$fam" "${MARKER[$fam]}" "$N" "$NR" "$NF" "$ALEN" >> "$OUT/03E_tree_pipeline_manifest.tsv"
done

echo "Fig3 03E four-family phylogeny pipeline PASS"
