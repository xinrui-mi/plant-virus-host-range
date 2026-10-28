# Checksums

`SHA256SUMS.txt` lists a SHA-256 for every file this release ships. Verify with:

```bash
cd <the released directory>
sha256sum -c checksums/SHA256SUMS.txt
```

The file is generated, never hand-edited:

```bash
bash checksums/generate_sha256sums.sh
```

## What is in scope

Every file **git tracks**, which is exactly the set a cloner receives. That deliberately excludes:

| Excluded | Why |
|---|---|
| `analysis/*/tables/`, `analysis/*/qc/` | Run artefacts. `04D`, for instance, writes 19 tables there; the released form of those numbers is `supplementary/Table_S3_Arthropod_vector/S3I..S3L`. |
| `analysis/*/pipeline_results/03L_FINAL_12family_host_spectrum/` | Output of `03L`; reproducible from the shipped trees in about 10 minutes. |
| `analysis/*/work/`, `analysis/*/plot/*/main|panels|supplementary/` | Scratch and rendered figure output. |
| `*.log`, `*.contree`, `*.iqtree`, `*.mldist`, `*.bionj`, `*.ckp.gz` | IQ-TREE side files. The published `.treefile` in `data/` is the artefact that matters. |

If you want to confirm that a regenerated output matches what the release describes, compare against
the `supplementary/` table rather than against a checksum: the checksum tells you the shipped file is
the one we shipped, not that your re-run agrees with it.

## Symlinked entries

Four entries in `data/analysis_data/` are **symlinks to files shipped elsewhere in the same release**,
not separate copies:

| Link | Target |
|---|---|
| `13_arthropod_sequence_supported_plant_ids.tsv.gz` | `arthropod_sequence_supported_union.tsv.gz` |
| `20_arthropod_sequence_evidence_FINAL.tsv.gz` | `arthropod_sequence_evidence_FINAL.tsv.gz` |
| `20_three_source_association_master_FINAL.tsv.gz` | `../analysis_input/association_analysis_input.tsv.gz` |
| `20_three_source_sequence_master_FINAL.tsv.gz` | `../analysis_input/sequence_analysis_input.tsv.gz` |

The aliases exist because the analysis scripts reference these three-source master files by their
historical cluster names, while the release ships the same bytes under descriptive names. Both paths
resolve inside the released tree, so no file is missing. `sha256sum -c` follows the links and verifies
the bytes; `find -type f` does not count them, which is why an extracted tree has 4 fewer regular
files than the 634 checksum lines.

## Stability warning

The list is only valid for the exact file contents at the moment of generation. Any edit after that
point makes it stale. Re-run the generator rather than editing the list by hand — a hand-edited
checksum list is worse than no checksum list, because it looks authoritative.
