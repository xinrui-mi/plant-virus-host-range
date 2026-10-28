#!/usr/bin/env bash
set -euo pipefail
source "${ANALYSIS_ROOT}/00_config/00_paths.sh"
BASE="${ANALYSIS_ROOT}/06_TSWV_evolution"
IN="${BASE}/pipeline_results/02_exact_haplotypes"
OUT="${BASE}/pipeline_results/03_alignment_and_phylogeny"
mkdir -p "$OUT" "${BASE}/logs"
printf 'gene\thaplotype_n\tprotein_alignment_length\ttrimmed_codon_length\ttreefile\n' > "$OUT/07E_five_gene_tree_manifest.tsv"
for gene in RdRp GnGc NSm N NSs; do
  D="$OUT/$gene"; mkdir -p "$D"
  FAA="$IN/${gene}.haplotypes.faa"; FNA="$IN/${gene}.haplotypes.fna"
  NSEQ=$(grep -c '^>' "$FAA" || true); [[ "$NSEQ" -ge 4 ]] || { echo "$gene has <4 exact haplotypes" >&2; exit 21; }
  if [[ "$NSEQ" -le 200 ]]; then
    mafft --globalpair --maxiterate 1000 --thread "$NCORES" "$FAA" > "$D/${gene}.protein.mafft.full.faa" 2> "$BASE/logs/${gene}.mafft.log"
  else
    mafft --auto --thread "$NCORES" "$FAA" > "$D/${gene}.protein.mafft.full.faa" 2> "$BASE/logs/${gene}.mafft.log"
  fi
  # Protein trimming and codon-aware backtranslation use the same trimAl rule.
  trimal -in "$D/${gene}.protein.mafft.full.faa" -out "$D/${gene}.protein.trim.faa" -automated1 > "$BASE/logs/${gene}.trimal.protein.log" 2>&1
  trimal -in "$D/${gene}.protein.mafft.full.faa" -out "$D/${gene}.codon.trim.fna" -automated1 -backtrans "$FNA" -ignorestopcodon > "$BASE/logs/${gene}.trimal.backtrans.log" 2>&1
  ALEN=$(awk 'BEGIN{s=""}/^>/{if(length(s)>0){print length(s);exit}next}{s=s$0}END{if(length(s)>0)print length(s)}' "$D/${gene}.protein.trim.faa" | head -1)
  NLEN=$(awk 'BEGIN{s=""}/^>/{if(length(s)>0){print length(s);exit}next}{s=s$0}END{if(length(s)>0)print length(s)}' "$D/${gene}.codon.trim.fna" | head -1)
  [[ $((NLEN % 3)) -eq 0 ]] || { echo "$gene trimmed codon alignment length not divisible by 3: $NLEN" >&2; exit 22; }
  cat > "$D/${gene}.codon_positions.partition.txt" <<EOF
DNA, codon1 = 1-${NLEN}\\3
DNA, codon2 = 2-${NLEN}\\3
DNA, codon3 = 3-${NLEN}\\3
EOF
  iqtree3 -s "$D/${gene}.codon.trim.fna" -p "$D/${gene}.codon_positions.partition.txt" -m MFP+MERGE -B 1000 --alrt 1000 -T "$NCORES" --seed 20260813 --prefix "$D/${gene}.codon_partitioned" > "$BASE/logs/${gene}.iqtree.stdout.log" 2> "$BASE/logs/${gene}.iqtree.stderr.log"
  printf '%s\t%s\t%s\t%s\t%s\n' "$gene" "$NSEQ" "$ALEN" "$NLEN" "$D/${gene}.codon_partitioned.treefile" >> "$OUT/07E_five_gene_tree_manifest.tsv"
done
echo "Fig7 five-gene alignment/tree pipeline PASS"
