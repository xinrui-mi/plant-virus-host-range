#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, os, re
from pathlib import Path
from Bio.Seq import Seq

REPO = Path(os.environ.get("REPO_ROOT", Path(os.environ["ANALYSIS_ROOT"]).parent))
ROOT = Path(os.environ["ANALYSIS_ROOT"])
SEQ_FINAL = Path(os.environ.get("SEQ_FINAL", REPO / "data/analysis_input/sequence_analysis_input.tsv.gz"))
ALL_FASTA = Path(os.environ.get("ALL_THREE_FASTA", REPO / "external/all_three_source.fasta"))
WORK = ROOT / "02_host_relatedness_and_phylogeny" / "work" / "marker_inputs"
WORK.mkdir(parents=True, exist_ok=True)
MIN_ORF_AA = int(os.environ.get("FIG3_RNA_MIN_ORF_AA", "200"))
FAMS = {"Potyviridae", "Tombusviridae"}


def base_acc(x): return re.sub(r"\.[0-9]+$", "", (x or "").strip())

def eligible(r):
    return r.get("virus_ictv_match_rank", "").strip().lower() in {"species", "genus", "family"} and r.get("virus_family", "") in FAMS

targets = {f: [] for f in FAMS}; fam_by_base = {}
_opener = gzip.open if str(SEQ_FINAL).endswith(".gz") else open
with _opener(SEQ_FINAL, "rt", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if eligible(r):
            sid = r["sequence_id"].strip(); fam = r["virus_family"]
            targets[fam].append(sid); fam_by_base[base_acc(sid)] = (sid, fam)


def fasta_iter(path):
    h=None; buf=[]
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                if h is not None: yield h, "".join(buf)
                h=line[1:].split()[0]; buf=[]
            else: buf.append(line.strip())
        if h is not None: yield h, "".join(buf)


def six_frame_orfs(nt):
    nt = re.sub(r"[^ACGTUNacgtun]", "N", nt).upper().replace("U", "T")
    seq = Seq(nt)
    out=[]
    for strand, s in [(1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            usable = ((len(s)-frame)//3)*3
            if usable < 3: continue
            aa = str(s[frame:frame+usable].translate(to_stop=False))
            cursor=0
            for part in aa.split("*"):
                if len(part) >= MIN_ORF_AA:
                    out.append((strand, frame, cursor, part))
                cursor += len(part)+1
    return out

handles = {fam: (WORK / f"{fam}.all_ORFs.faa").open("w") for fam in FAMS}
maps = {fam: (WORK / f"{fam}.all_ORFs_map.tsv").open("w", newline="") for fam in FAMS}
writers={}
for fam in FAMS:
    writers[fam]=csv.writer(maps[fam], delimiter="\t")
    writers[fam].writerow(["orf_id","sequence_id","strand","frame","aa_length"])
found={fam:set() for fam in FAMS}; orf_n={fam:0 for fam in FAMS}
try:
    for hid, nt in fasta_iter(ALL_FASTA):
        key=base_acc(hid)
        if key not in fam_by_base: continue
        sid,fam=fam_by_base[key]; found[fam].add(sid)
        for idx,(strand,frame,offset,aa) in enumerate(six_frame_orfs(nt),1):
            oid=f"{sid}__ORF{idx:04d}"
            handles[fam].write(f">{oid}\n{aa}\n")
            writers[fam].writerow([oid,sid,strand,frame,len(aa)]);orf_n[fam]+=1
finally:
    for x in handles.values(): x.close()
    for x in maps.values(): x.close()

with (WORK / "03E_RNA_ORF_QC.tsv").open("w", newline="") as f:
    w=csv.writer(f, delimiter="\t");w.writerow(["family","target_sequence_n","sequence_recovered_n","missing_sequence_n","orf_n","min_orf_aa"])
    for fam in sorted(FAMS):
        missing=set(targets[fam])-found[fam]
        w.writerow([fam,len(targets[fam]),len(found[fam]),len(missing),orf_n[fam],MIN_ORF_AA])
        (WORK/f"{fam}.missing_from_all_three_fasta.txt").write_text("\n".join(sorted(missing))+("\n" if missing else ""))
        print(fam,"targets",len(targets[fam]),"recovered",len(found[fam]),"ORFs",orf_n[fam],"missing",len(missing))
