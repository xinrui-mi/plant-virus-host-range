#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Parse CD-HIT .clstr into deterministic membership TSV."""

from __future__ import annotations
import argparse, csv, re
from pathlib import Path

cluster_re = re.compile(r"^>Cluster\s+(\d+)")
member_re = re.compile(
    r"^\s*\d+\s+\d+aa,\s+>(.+?)\.\.\.\s+(.*)$"
)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--clstr", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    family = args.family
    clstr = Path(args.clstr)
    out = Path(args.out)

    rows = []
    current = None
    cluster_members = {}
    with clstr.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            m = cluster_re.match(line)
            if m:
                current = int(m.group(1))
                cluster_members[current] = []
                continue
            m = member_re.match(line)
            if m and current is not None:
                sid = m.group(1)
                tail = m.group(2).strip()
                is_rep = tail == "*"
                identity = "100" if is_rep else ""
                mi = re.search(r"at\s+([0-9.]+)%", tail)
                if mi:
                    identity = mi.group(1)
                cluster_members[current].append((sid, is_rep, identity))

    for cid, members in sorted(cluster_members.items()):
        reps = [x[0] for x in members if x[1]]
        if len(reps) != 1:
            raise SystemExit(f"{family} cluster {cid} has {len(reps)} representatives")
        rep = reps[0]
        cluster_id = f"{family}_AA95_{cid+1:05d}"
        for sid, is_rep, identity in members:
            rows.append({
                "virus_family": family,
                "cluster95_id": cluster_id,
                "representative_sequence_id": rep,
                "member_sequence_id": sid,
                "is_representative": "TRUE" if is_rep else "FALSE",
                "identity_to_representative_pct": identity,
            })

    out.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "virus_family", "cluster95_id", "representative_sequence_id",
        "member_sequence_id", "is_representative",
        "identity_to_representative_pct"
    ]
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(
        f"CD-HIT membership PASS {family}: "
        f"clusters={len(cluster_members)} members={len(rows)}"
    )

if __name__ == "__main__":
    main()
