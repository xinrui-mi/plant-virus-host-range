#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, re

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("clstr")
    ap.add_argument("output")
    args = ap.parse_args()

    rows = []
    cluster = None
    members = []
    rep = None

    def flush():
        nonlocal rows, cluster, members, rep
        if cluster is None:
            return
        if rep is None:
            raise RuntimeError(f"No representative in {cluster}")
        for i, m in enumerate(members, start=1):
            rows.append({
                "cdhit95_cluster_id": cluster,
                "representative_sequence_id": rep,
                "member_sequence_id": m,
                "member_index": i,
                "is_representative": "YES" if m == rep else "NO",
                "cluster_member_n": len(members),
            })

    with open(args.clstr) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">Cluster "):
                flush()
                cluster = "AA95_" + line.split()[-1]
                members = []
                rep = None
                continue

            m = re.search(r">(.+?)\.\.\.", line)
            if not m:
                continue
            sid = m.group(1)
            members.append(sid)
            if line.endswith("*"):
                rep = sid

    flush()

    fields = [
        "cdhit95_cluster_id","representative_sequence_id",
        "member_sequence_id","member_index",
        "is_representative","cluster_member_n"
    ]
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"cluster_membership_rows={len(rows)}")
    print(f"cluster_n={len(set(r['cdhit95_cluster_id'] for r in rows))}")

if __name__ == "__main__":
    main()
