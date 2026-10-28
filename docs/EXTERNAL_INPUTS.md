# External inputs you must obtain yourself

The released reproducibility chain (levels L2 and L3 in
`docs/reproducibility.md`) needs **none** of what follows. It starts from the
released analysis-ready datasets and the frozen derived inputs shipped in
`data/`.

This page exists for two other purposes:

1. you want to **rebuild** the family phylogenies or the TSWV gene-level
   analyses from raw sequences;
2. you want to understand exactly where every released number came from, so you
   can audit the upstream provenance.

---

## 1. Source databases (needed only to rebuild from scratch)

| Resource | Role in the study | How to obtain |
|---|---|---|
| NCBI Virus | primary virus sequence + host-association source | NCBI Virus portal / E-utilities |
| Virus-Host DB | secondary host-association source | <https://www.genome.jp/virushostdb/> (CC BY-NC-SA 4.0 — see §5) |
| Land Flora RNA virome | plant RNA-virome records | the source study's released dataset |
| VIP-DB | documented arthropod vector relationships; transmission modes | the VIP-DB publication and its accompanying table |
| ICTV | current taxonomy (realm to species), used to define `virus_family` / `virus_genus` eligibility | ictv.global |
| EPPO / FAO-IPPC | regulatory records used in the manually curated external-evidence layer | eppo.int; ippc.int |

The harmonisation that turns the raw downloads into the released tables is
**outside the repository reproducibility scope**: it was interactive, it spans a
local R session and an HPC cluster, and it embeds manual curation decisions. It
is not reconstructed here, and no substitute implementation is provided, because
a re-implementation would yield a *different* database rather than a check on
the published one.

---

## 2. Third-party files referenced by the optional pipeline scripts

Set `EXTERNAL_DATA_DIR` and place these inside it:

```bash
export EXTERNAL_DATA_DIR=/path/to/external
source analysis/00_config/00_paths.sh
```

| Variable | File | Expected location | Used by |
|---|---|---|---|
| `PFAM_HMM` | `Pfam-A.hmm` | `$EXTERNAL_DATA_DIR/Pfam-A.hmm` | `03E_extract_dna_markers_fullpfam.py`, `03K_A_hmmscan_one_family_FIXED.sh` |
| `DNA_GENBANK` | `dna_allaccession.gb` | `$EXTERNAL_DATA_DIR/dna_allaccession.gb` | `03E_extract_dna_cds.py` |
| `ALL_THREE_FASTA` | `all_three_source.fasta` | `$EXTERNAL_DATA_DIR/all_three_source.fasta` | `03D_prepare_marker_inputs.py`, `03H_*`, `03J_*`, `03K_A_*` |
| `RDRP_CORE_FASTA` | `08_RdRP_core.faa` | `$EXTERNAL_DATA_DIR/08_RdRP_core.faa` | `03H_*`, `03I_prepare_RNA_family_RdRP_core.py` |
| `TSWV_SOURCE_GB` | `sequence.gb` | `$EXTERNAL_DATA_DIR/TSWV_final_pipeline/input/sequence.gb` | `07B_TSWV_membership_audit.py`, `07C_extract_TSWV_genes.py` |
| `VIPDB_RAW` | `VIPDB_relation_raw.xlsx` | `$EXTERNAL_DATA_DIR/VIPDB_relation_raw.xlsx` | `04A_VIPDB_normalize.R`, `additional_validation/transmission_mode_consistency/run_specificity_control.sh` |

### Digests of the files used in the frozen run

Captured at freeze time (`analysis/00_config/input_lock/00_input_SHA256.tsv`).
Compare with `sha256sum <file>` before rebuilding: HMMER searches and flat-file
parsing are version sensitive, and a different HMM release will not reproduce
the marker sets.

| File | SHA256 |
|---|---|
| `Pfam-A.hmm` | `a78a42d6faf265b6bfca59e8f062d06fae6083ce2c6e335d7b381f20b82b7903` |
| `dna_allaccession.gb` | `1430ce5e69d5ca7d74dbcd34b86295d069b7b162edcdf9ef868303293b48aa23` |
| `sequence.gb` (TSWV) | `bc8e6f10696b26dd99312bea9291f84d0bdec4cb7db375293e9fe12999414ff5` |
| `VIPDB_relation_raw.xlsx` | `77ce475c544d5b8a9b5f133a395ddc132cccffb16b73bcc68b81679cc1e20e22` |
| `08_RdRP_core.faa` | `02a32623fe469749c43624f8fcc4918a754c0403614ea246ced5748c05062a71` |
| `all_three_source.fasta` | `dfd0d78ac1453b301e3bc3f9f785bf8fecb004e88c3a6d8ebaf4ceaca055daee` |
| final `20_three_source_sequence_master_FINAL.tsv` | `0e6224455dc3c0983334dd6bf2bedd870c4169b876eb2e4596fcdaaa8929a0f6` |
| final `20_three_source_association_master_FINAL.tsv` | `d9d40aff9692ed9c2f4dfdc91c676d64d9cca3c2f134f4b99e6de1f0d13b7620` |
| final `20_arthropod_sequence_evidence_FINAL.tsv` | `93c2b2096272bc6eb060dea3115fb8f109a635f9bea61712c5f09913124b176f` |
| `candidate_external_evidence_legacy_traceable.tsv` | `78180e1aa60e2affd6c02861cd7a8f932b518abf3db049c9c948f8545c649164` |

The first six are **third-party** files, not redistributed here. The last four
are project-produced; the analysis-ready projections of the two master tables
*are* redistributed, as `data/analysis_input/*.tsv.gz`.

---

## 3. Required software for the optional pipelines only

`cd-hit`, `mafft`, `trimal`, `iqtree3`, `hmmsearch`, `hmmfetch`, `hyphy`
(plus `HYPHYMPI` / `mpirun` for MPI runs). Versions used are in
`environment/software_versions.txt`.

`analysis/00_config/00_preflight.sh` reports whether they are present and does
**not** fail when they are missing, precisely because they are optional.

---

## 4. What happens when a third-party file is unavailable

The released code degrades explicitly rather than silently:

| Script | Behaviour without its third-party input |
|---|---|
| `04A_VIPDB_normalize.R` | Restores the released normalized VIPDB tables into `<module>/tables/` from `data/analysis_data/` and exits successfully. Downstream modules then run unchanged. |
| `07B_TSWV_membership_audit.py` | Reports the TSWV membership found in the released sequence input, writes the GATE file with reason `TSWV_SOURCE_GB_NOT_AVAILABLE`, and stops. |
| `03E_*`, `03H_*`, `03I_*`, `03J_*`, `03K_*`, `07C`-`07M` | Not part of the released chain; their frozen outputs are shipped instead. |
| `additional_validation/.../04D_transmission_specificity_control.R` | Requires `VIPDB_RAW`; no fallback. This module is `NOT_IN_PAPER`, so nothing in the manuscript depends on it. |

---

## 5. Licensing and redistribution — read before reuse

*(Status of this section: the statements below describe what the repository
does. Upstream terms marked **VERIFIED** below were checked against the
upstream source on 2026-10-27, with the verification URL given in each entry.
No term is asserted beyond what the upstream source states.)*

| Layer | Licence | Notes |
|---|---|---|
| Repository code (`analysis/`, `additional_validation/`, scripts) | MIT — see `LICENSE` | |
| Project-produced derived tables and compiled datasets (`data/`, `supplementary/`, `docs/`) | CC BY 4.0 — see `LICENSE_DATA` | Applies to the compilation, selection and derived-value layer produced by this project. |
| Manual curation content (`data/external_evidence/`) | CC BY 4.0 | The curation and provenance annotations are project output; the underlying facts belong to EPPO / FAO-IPPC / VIP-DB / the cited literature. |
| NCBI Virus records | **VERIFIED — redistributable under NCBI policy** | Checked 2026-10-27: <https://www.ncbi.nlm.nih.gov/home/about/policies/> states "NCBI itself places no restrictions on the use or distribution of the data contained therein". Caveat kept: NCBI does not evaluate or waive any third-party claims (e.g. patents or copyright) that individual submitters may assert on specific sequences. Redistributing *derived, analysis-ready* projections is not re-licensing the records; `LICENSE_DATA` does not re-licence them. |
| Virus-Host DB records | **VERIFIED — CC BY-NC-SA 4.0 (incompatible with CC BY 4.0)** | Checked 2026-10-27: <https://www.genome.jp/virushostdb/about> states the database is licensed CC BY-NC-SA 4.0. Redistribution of records is permitted with attribution, for **non-commercial** use, share-alike. **The NC (and SA) terms are not compatible with this package's CC BY 4.0**: the embedded Virus-Host DB records in the released tables remain under CC BY-NC-SA 4.0 and are expressly *excluded* from the package's CC BY 4.0 grant. |
| Land Flora RNA virome records | **VERIFIED — open publication; redistributable with citation** | Checked 2026-10-27: this corresponds to **Liu et al. 2026, "The Origins and Evolution of the RNA Virome in Land Flora", *Molecular Plant* 19(10), doi:10.1016/j.molp.2026.07.023** (Liu, Gong, Koonin, Dolja & Han). Its data-availability statement releases all raw data, intermediate files and final result tables on figshare: **doi:10.6084/m9.figshare.30688322** (analysis scripts: github.com/liuqiyan67/Plant_virome). The 2,516 records in the curated table carry only this project's internal identifiers (`TS_nucl_*`); they are derived from that openly published dataset. Redistribute with citation of the paper and the figshare DOI; confirm the licence shown on the figshare item page before any further redistribution. |
| VIP-DB documented vector table | **VERIFIED — open publication; normalized projection redistributable with citation** | Checked 2026-10-27: VIP-DB (Virus–Insect–Plant Database), Deng et al. 2026, *Viruses* 18(6):679, doi:10.3390/v18060679 (MDPI, open access; standard MDPI CC BY terms apply to the paper/dataset); the resource is publicly downloadable from <http://viplant.cn/> and reports **1,375 integrated virus–insect–plant records**, exactly matching the 1,375 rows of the released `VIPDB_relation_normalized.tsv.gz` — i.e. the shipped projection is the published open dataset. The raw workbook itself is still **not** redistributed here; obtain it from the source and cite the paper. |
| ICTV taxonomy | **VERIFIED — CC BY 4.0 (compatible with this package)** | Checked 2026-10-27: <https://ictv.global> states "Unless otherwise noted, this work is licensed under the CC BY 4.0"; the MSL/VMR releases carry Zenodo DOIs under the same licence. Compatible with `LICENSE_DATA`. |
| Pfam-A HMM | **VERIFIED — CC0 1.0; not redistributed here** | Checked 2026-10-27: <https://www.ebi.ac.uk/interpro/about/license/> — InterPro/Pfam downloadable data are released under CC0 1.0 (public domain). Obtain directly from EBI; pin the digest below. |
| GenBank flat files | **VERIFIED — same terms as NCBI Virus records above** | Subject to NCBI policy plus any submitter claims; not redistributed here. |

> **Plain statement for readers.** `LICENSE_DATA` (CC BY 4.0) covers what this
> project produced. It does **not** re-licence any third-party source data, and
> nothing in this repository should be read as claiming that NCBI Virus,
> Virus-Host DB, Land Flora RNA virome or VIP-DB content is available under
> CC BY 4.0. In particular, the embedded **Virus-Host DB records remain under
> CC BY-NC-SA 4.0**, which is *not* compatible with CC BY 4.0 — reuse them
> under their own terms. Upstream reuse conditions above were verified on
> 2026-10-27 against the listed sources; if you redistribute a subset of the
> derived tables, comply with the upstream terms above and cite the upstream
> sources.

### 5a. Three tiers of content in this repository

| Tier | What it is | Licence status |
|---|---|---|
| **A. Author-generated statistics and code** | every script; the fitted coefficients, p-values, q-values, permutation and bootstrap distributions in `supplementary/` and `analysis/*/tables/`; the curated evidence reading in `data/external_evidence/`; the compilation/selection layer of the curated datasets | CC BY 4.0 (`LICENSE_DATA`) / MIT (`LICENSE`) — this is the project's own work |
| **B. Derived tables that embed third-party records** | `data/curated/*.tsv.gz` (native fields: TaxIDs, taxonomy, host names from NCBI Virus / Virus-Host DB / Land Flora); `data/analysis_data/VIPDB_*.tsv.gz` (normalized VIP-DB projection); ICTV-derived taxonomy columns throughout | The *derivation* is CC BY 4.0; the embedded third-party records remain under tier C terms. `source` / `source_record_id` columns let you separate the two. |
| **C. Third-party raw sources** | NCBI Virus / GenBank records, Virus-Host DB, the Land Flora RNA virome study dataset, the VIP-DB raw workbook, ICTV taxonomy, EPPO and FAO/IPPC documents, Pfam-A HMM | NOT redistributed as raw dumps; each remains under its own (now verified — see the table above) terms |

### 5b. Redistribution status of the shipped files that embed third-party records

These files are shipped because the analysis needs them. The upstream terms for
each embedded source were verified on 2026-10-27 (see the table in §5 for the
verification URLs); the status column states what that means for
redistribution.

| File(s) in this release | Third-party content embedded | Upstream | Status |
|---|---|---|---|
| `data/curated/plant_virus_sequence_dataset.tsv.gz` | accessions, TaxIDs, taxonomy, lengths | NCBI Virus / GenBank | OK to redistribute — NCBI places no restrictions on its data (checked 2026-10-27); individual submitter claims, if any, follow the sequences. Cite NCBI/GenBank accessions. |
| `data/curated/virus_host_association_dataset.tsv.gz` | association records + host names | NCBI Virus; Virus-Host DB | OK to redistribute with two parallel terms: NCBI-derived content unrestricted; the **Virus-Host DB records remain CC BY-NC-SA 4.0** (non-commercial, share-alike) and are excluded from the package's CC BY 4.0 grant. Cite both sources. |
| `data/analysis_data/VIPDB_relation_normalized.tsv.gz`, `VIPDB_virus_level_summary.tsv.gz` | normalized VIP-DB vector relations (1,375 records) | VIP-DB (Deng et al. 2026, *Viruses* 18:679; viplant.cn) | OK to redistribute — the projection matches the openly published dataset; cite the paper and viplant.cn. |
| taxonomy columns (`virus_realm`…`virus_species`) in all curated/derived tables | ICTV taxonomy | ICTV | OK to redistribute — ICTV data are CC BY 4.0 (checked 2026-10-27), compatible with `LICENSE_DATA`. Cite ICTV. |
| `data/analysis_data/` arthropod evidence tables | build on the above three sources | inherited | Licence status is inherited from the now-confirmed upstream terms above: the derivation layer is CC BY 4.0, embedded records keep their upstream terms (NCBI unrestricted; Virus-Host DB CC BY-NC-SA 4.0). |
| Land Flora RNA virome records (inside the curated sequence table, `source` = Land Flora) | 2,516 sequence records | Liu et al. 2026, *Molecular Plant* 19(10), doi:10.1016/j.molp.2026.07.023; data on figshare doi:10.6084/m9.figshare.30688322 | OK to redistribute — records derive from the authors' openly published figshare dataset (checked 2026-10-27). Cite the paper and the figshare DOI. |
| — not shipped — | Pfam-A HMM (CC0 1.0); GenBank flat files (NCBI policy); VIP-DB raw workbook (open publication); TSWV source GenBank (NCBI policy) | Pfam, NCBI, VIP-DB | excluded from redistribution; obtain from the source and pin the digests in §2 |

---

## 6. Cite these too

Republishing derived data obliges you to acknowledge the upstream resources:

NCBI Virus · Virus-Host DB · Land Flora RNA virome · VIP-DB · ICTV ·
EPPO / FAO-IPPC (records referenced in `data/external_evidence/`) ·
Pfam · GenBank.
