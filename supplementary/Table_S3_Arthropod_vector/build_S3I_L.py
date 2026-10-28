#!/usr/bin/env python3
"""Rebuild Table S3 sections I-L from the 04D output tables.

`04D_transmission_biology_extended_v3.R` writes 19 result tables into
`analysis/03_arthropod_transmission_ecology/tables/` (20_04D_* ... 38_04D_*).
Table S3 sections I-L are stacked views of those tables: each section keeps
every row of its source tables and adds two leading columns that say which
source table a row came from and which block it belongs to, so a reader can
trace a number in the supplementary table straight back to a table name.

Because the blocks have different column sets (9 to 41 columns), each block is
written as a self-contained sub-table with its own header row rather than
merged into one wide rectangle. Column meaning is therefore never ambiguous and
no cell is ever padded with a placeholder.

Nothing is recomputed and no value is edited: every data cell is copied
verbatim from the 04D output.
"""

import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(REPO, "analysis", "03_arthropod_transmission_ecology", "tables")
DST = os.path.join(REPO, "supplementary", "Table_S3_Arthropod_vector")

# section -> (output file name, [(source table, block label)])
SECTIONS = {
    "S3I_vector_family_by_virus_genus.tsv": [
        ("24_04D_vector_family_by_virus_genus.tsv", "vector_family_x_virus_genus"),
    ],
    "S3J_five_group_host_range_models.tsv": [
        ("20_04D_five_group_descriptive_summary.tsv", "descriptive_summary"),
        ("25_04D_five_group_host_range_model_terms.tsv", "host_range_model_terms"),
        ("26_04D_five_group_host_range_primary_contrasts.tsv", "host_range_primary_contrasts"),
        ("27_04D_lineage_leave_one_out_primary_contrasts.tsv", "lineage_leave_one_out"),
        ("34_04D_lineage_driver_ranked_summary.tsv", "lineage_driver_ranking"),
    ],
    # S3K = within-mode lineage leave-out sensitivity (matches Data_S1.xlsx
    # section K1-K3); S3L = country breadth and evidence models (Excel L1-L3).
    # The previous release had these two file names swapped relative to the
    # Excel/Manuscript labelling; fixed 2026-10-27 (UNRESOLVED_ITEMS.md U25).
    "S3K_within_mode_lineage_sensitivity.tsv": [
        ("23_04D_within_family_mode_overlap_audit.tsv", "within_family_mode_overlap_audit"),
        ("35_04D_leave_one_lineage_out_eligibility_audit.tsv", "leave_one_lineage_out_eligibility"),
        ("36_04D_within_mode_lineage_leaveout_primary_contrasts.tsv", "within_mode_leaveout_contrasts"),
        ("37_04D_within_mode_lineage_influence_ranked_summary.tsv", "within_mode_influence_ranking"),
        ("38_04D_within_mode_lineage_eligibility_audit.tsv", "within_mode_eligibility"),
    ],
    "S3L_country_breadth_and_evidence_models.tsv": [
        ("28_04D_country_breadth_model_terms.tsv", "country_breadth_model_terms"),
        ("29_04D_country_breadth_primary_contrasts.tsv", "country_breadth_primary_contrasts"),
        ("30_04D_S0S2_vs_VIPDB_raw_fisher.tsv", "S0S2_vs_VIPDB_raw_fisher"),
        ("31_04D_S0S2_vs_VIPDB_adjusted_model_terms.tsv", "S0S2_vs_VIPDB_adjusted_terms"),
        ("32_04D_S0S2_vs_VIPDB_adjusted_primary.tsv", "S0S2_vs_VIPDB_adjusted_primary"),
        ("33_04D_binary_evidence_country_breadth_primary.tsv", "binary_evidence_country_breadth"),
    ],
}

# Per-table column counts differ (41 columns for the leave-one-out blocks, 9-18
# for the rest). Each table block is written as its own contiguous run so that
# column meaning is never ambiguous; rows are padded to the widest block in the
# section and the padding is left empty rather than filled with a placeholder.


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        text = fh.read()
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    header = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:]]
    return header, rows


def main():
    missing = []
    for out_name, blocks in SECTIONS.items():
        loaded = []
        for fname, label in blocks:
            p = os.path.join(SRC, fname)
            if not os.path.exists(p):
                missing.append(fname)
                continue
            header, rows = read_tsv(p)
            loaded.append((fname, label, header, rows))
        if missing:
            continue

        wide = max(len(h) for _, _, h, _ in loaded)
        out_path = os.path.join(DST, out_name)
        n_data = 0
        with open(out_path, "w", encoding="utf-8", newline="") as fh:
            for fname, label, header, rows in loaded:
                # One header row per block: the blocks genuinely have different
                # column sets, so a single shared header would be wrong.
                top = ["section_block", "source_table"] + header
                top += [""] * (wide + 2 - len(top))
                fh.write("\t".join(top) + "\n")
                for r in rows:
                    line = [label, fname] + r
                    line += [""] * (wide + 2 - len(line))
                    fh.write("\t".join(line) + "\n")
                    n_data += 1
        print(f"{out_name}: {len(loaded)} block(s), {n_data} data rows, {wide + 2} columns")


    if missing:
        print("MISSING SOURCE TABLES:", ", ".join(sorted(set(missing))), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
