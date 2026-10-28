#!/usr/bin/env python3
from pathlib import Path
import argparse
import csv

STD_TABLE = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

# Translation-table-4-like audit for Mitoviridae.
MITO4_TABLE = dict(STD_TABLE)
MITO4_TABLE['TGA'] = 'W'
MITO4_TABLE['ATA'] = 'M'
MITO4_TABLE['AGA'] = 'S'
MITO4_TABLE['AGG'] = 'S'

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(s):
    return s.translate(COMP)[::-1]

def read_fasta(path):
    h = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(seq).upper().replace("U","T")
                h = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if h is not None:
        yield h, "".join(seq).upper().replace("U","T")

def translate_frame(nt, frame, table):
    aa = []
    for i in range(frame, len(nt)-2, 3):
        aa.append(table.get(nt[i:i+3], "X"))
    return "".join(aa)

def stop_to_stop_orfs(nt, strand, table_name, table, min_aa):
    seq = nt if strand == "+" else revcomp(nt)
    out = []

    for frame in range(3):
        aa = translate_frame(seq, frame, table)
        start = 0
        for j in range(len(aa)+1):
            end_here = (j == len(aa)) or (j < len(aa) and aa[j] == "*")
            if not end_here:
                continue

            frag = aa[start:j]
            if len(frag) >= min_aa:
                nt_start_local = frame + start * 3
                nt_end_local = frame + j * 3

                if strand == "+":
                    nt_start = nt_start_local + 1
                    nt_end = nt_end_local
                else:
                    nt_start = len(nt) - nt_end_local + 1
                    nt_end = len(nt) - nt_start_local

                out.append({
                    "strand": strand,
                    "frame": frame + 1,
                    "translation_table": table_name,
                    "nt_start": min(nt_start, nt_end),
                    "nt_end": max(nt_start, nt_end),
                    "aa": frag,
                    "aa_length": len(frag),
                })
            start = j + 1

    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--out-faa", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--min-aa", type=int, default=120)
    ap.add_argument("--max-orfs-per-seq", type=int, default=12)
    args = ap.parse_args()

    fasta_rows = []
    rows = []
    seq_total = 0
    seq_with_orf = 0

    for sid, nt in read_fasta(args.input):
        seq_total += 1
        cand = []
        cand.extend(stop_to_stop_orfs(nt, "+", "standard", STD_TABLE, args.min_aa))
        cand.extend(stop_to_stop_orfs(nt, "-", "standard", STD_TABLE, args.min_aa))

        if args.family == "Mitoviridae":
            cand.extend(stop_to_stop_orfs(
                nt, "+", "mitochondrial_code4_audit", MITO4_TABLE, args.min_aa
            ))
            cand.extend(stop_to_stop_orfs(
                nt, "-", "mitochondrial_code4_audit", MITO4_TABLE, args.min_aa
            ))

        cand.sort(key=lambda x: x["aa_length"], reverse=True)
        cand = cand[:args.max_orfs_per_seq]

        if cand:
            seq_with_orf += 1

        for rank, x in enumerate(cand, start=1):
            oid = (
                f"{sid}|orf{rank}|{x['translation_table']}|"
                f"{x['strand']}{x['frame']}|{x['nt_start']}-{x['nt_end']}"
            )
            fasta_rows.append((oid, x["aa"]))
            rows.append({
                "virus_family": args.family,
                "sequence_id": sid,
                "orf_id": oid,
                "orf_rank_by_length": rank,
                "translation_table": x["translation_table"],
                "strand": x["strand"],
                "frame": x["frame"],
                "nt_start": x["nt_start"],
                "nt_end": x["nt_end"],
                "aa_length": x["aa_length"],
            })

    Path(args.out_faa).parent.mkdir(parents=True, exist_ok=True)

    with open(args.out_faa, "w") as f:
        for h, s in fasta_rows:
            f.write(f">{h}\n{s}\n")

    fields = [
        "virus_family","sequence_id","orf_id","orf_rank_by_length",
        "translation_table","strand","frame","nt_start","nt_end","aa_length"
    ]
    with open(args.out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(
        f"{args.family}\tinput_sequence_n={seq_total}\t"
        f"sequence_with_candidate_orf_n={seq_with_orf}\t"
        f"candidate_orf_n={len(rows)}"
    )

if __name__ == "__main__":
    main()
