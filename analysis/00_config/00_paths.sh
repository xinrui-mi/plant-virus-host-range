#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Source this file from the repository root before running any analysis module:
#
#     source analysis/00_config/00_paths.sh
#
# It exports the roots and every input path used by the released scripts.
# Nothing in this file is machine specific: all defaults resolve inside the
# repository, so a fresh clone runs without editing any path.
#
# Upstream databases that cannot be redistributed (NCBI Virus dumps,
# Virus-Host DB dumps, the raw VIPDB workbook, Pfam HMM, source GenBank
# records) are referenced only through EXTERNAL_DATA_DIR and are needed
# solely by the tree-rebuilding pipelines.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Locate this file, and therefore the repository, without assuming any shell.
#
# The obvious `dirname "${BASH_SOURCE[0]}"` is wrong when the file is sourced
# from a shell that has no BASH_SOURCE (zsh, or `sh -c`): the expansion is empty,
# `dirname ""` is ".", and `${...}/..` resolves to the *parent* of the current
# directory. That yields a wrong ANALYSIS_ROOT — one level too high — which still
# exists, so nothing downstream fails loudly.
#
# So: try BASH_SOURCE, then $0, then a relative probe, and finally validate that
# the candidate actually looks like this repository. If none of them do, stop
# rather than export a plausible-looking wrong path.
# ---------------------------------------------------------------------------
_is_repo_analysis_root() {
  # A candidate is valid only if it contains this very file and the config dir.
  [ -f "$1/00_config/00_paths.sh" ] && [ -d "$1/00_config" ]
}

_ANALYSIS_ROOT=""
# ${BASH_SOURCE[0]:-} is a bash-ism; under dash/sh it is a syntax error, not an
# empty expansion. Split the two cases so the file also survives `. file` in a
# POSIX shell, where $0 is the only clue available.
if [ -n "${BASH_VERSION:-}" ]; then
  _src_candidates="${BASH_SOURCE[0]:-} $0"
else
  _src_candidates="$0"
fi
for _cand_src in $_src_candidates; do
  [ -n "$_cand_src" ] || continue
  [ -f "$_cand_src" ] || continue
  _cand="$(cd "$(dirname "$_cand_src")/.." 2>/dev/null && pwd)" || continue
  if _is_repo_analysis_root "$_cand"; then _ANALYSIS_ROOT="$_cand"; break; fi
done
unset _src_candidates
if [ -z "$_ANALYSIS_ROOT" ] && [ -f "analysis/00_config/00_paths.sh" ]; then
  _ANALYSIS_ROOT="$(cd "analysis" && pwd)"
fi
unset _cand_src _cand

if [ -z "$_ANALYSIS_ROOT" ] || ! _is_repo_analysis_root "$_ANALYSIS_ROOT"; then
  echo "00_paths.sh: cannot locate the repository analysis root." >&2
  echo "  Source this file as:  source <repo>/analysis/00_config/00_paths.sh" >&2
  echo "  or set ANALYSIS_ROOT explicitly before sourcing." >&2
  return 1 2>/dev/null || exit 1
fi

export ANALYSIS_ROOT="${ANALYSIS_ROOT:-${_ANALYSIS_ROOT}}"
export REPO_ROOT="${REPO_ROOT:-$(cd "${ANALYSIS_ROOT}/.." && pwd)}"

export DATA_DIR="${DATA_DIR:-${REPO_ROOT}/data}"
export ANALYSIS_DATA="${ANALYSIS_DATA:-${DATA_DIR}/analysis_data}"

# --- Released analysis-ready inputs -----------------------------------------
export SEQ_FINAL="${SEQ_FINAL:-${DATA_DIR}/analysis_input/sequence_analysis_input.tsv.gz}"
export ASSOC_FINAL="${ASSOC_FINAL:-${DATA_DIR}/analysis_input/association_analysis_input.tsv.gz}"
export ARTH_SEQ_FINAL="${ARTH_SEQ_FINAL:-${ANALYSIS_DATA}/arthropod_sequence_evidence_FINAL.tsv.gz}"

# --- Released derived analysis inputs ---------------------------------------
export VIPDB_NORMALIZED="${VIPDB_NORMALIZED:-${ANALYSIS_DATA}/VIPDB_relation_normalized.tsv.gz}"
export VIPDB_VIRUS_LEVEL="${VIPDB_VIRUS_LEVEL:-${ANALYSIS_DATA}/VIPDB_virus_level_summary.tsv.gz}"
export CANDIDATE_EVIDENCE="${CANDIDATE_EVIDENCE:-${DATA_DIR}/external_evidence/candidate_external_evidence.tsv.gz}"

export PHYLO_DIR="${PHYLO_DIR:-${ANALYSIS_DATA}/phylogeny}"
# 12 family-specific marker trees + AA95 membership and tip-to-host tables
export FAMILY_TREE_DIR="${FAMILY_TREE_DIR:-${PHYLO_DIR}/twelve_family}"
export TSWV_TREE_DIR="${TSWV_TREE_DIR:-${PHYLO_DIR}/TSWV}"
export TSWV_LMS_META="${TSWV_LMS_META:-${TSWV_TREE_DIR}/LMS_matched_isolate_metadata.tsv}"

# --- Non-redistributable upstream inputs (tree rebuilding only) -------------
# EXTERNAL_DATA_DIR stays EMPTY unless the user points it at a local copy of the
# upstream databases. The dependent variables are only expanded when it is set;
# leaving them empty (rather than "/file.xlsx" or "./file.xlsx") lets the
# pipelines detect "not supplied" cleanly and fall back to released products.
export EXTERNAL_DATA_DIR="${EXTERNAL_DATA_DIR:-}"
if [ -n "${EXTERNAL_DATA_DIR}" ]; then
  export ALL_THREE_FASTA="${ALL_THREE_FASTA:-${EXTERNAL_DATA_DIR}/all_three_source.fasta}"
  export RDRP_CORE_FASTA="${RDRP_CORE_FASTA:-${EXTERNAL_DATA_DIR}/08_RdRP_core.faa}"
  export VIPDB_RAW="${VIPDB_RAW:-${EXTERNAL_DATA_DIR}/VIPDB_relation_raw.xlsx}"
  export DNA_GENBANK="${DNA_GENBANK:-${EXTERNAL_DATA_DIR}/dna_allaccession.gb}"
  export PFAM_HMM="${PFAM_HMM:-${EXTERNAL_DATA_DIR}/Pfam-A.hmm}"
  export TSWV_SOURCE_GB="${TSWV_SOURCE_GB:-${EXTERNAL_DATA_DIR}/TSWV_final_pipeline/input/sequence.gb}"
else
  export ALL_THREE_FASTA="${ALL_THREE_FASTA:-}"
  export RDRP_CORE_FASTA="${RDRP_CORE_FASTA:-}"
  export VIPDB_RAW="${VIPDB_RAW:-}"
  export DNA_GENBANK="${DNA_GENBANK:-}"
  export PFAM_HMM="${PFAM_HMM:-}"
  export TSWV_SOURCE_GB="${TSWV_SOURCE_GB:-}"
fi

export NCORES="${NCORES:-12}"
