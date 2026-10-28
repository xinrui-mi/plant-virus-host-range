#!/usr/bin/env python3
from pathlib import Path
import argparse, re

VALID = set("ACDEFGHIKLMNPQRSTVWYX")

def read_fasta(path):
    h = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(buf)
                h = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if h is not None:
        yield h, "".join(buf)

def clean_seq(s):
    s = s.upper().replace("*", "")
    return "".join(ch if ch in VALID else "X" for ch in s)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output")
    args = ap.parse_args()

    seen = set()
    n = 0
    with open(args.output, "w") as out:
        for h, s in read_fasta(args.input):
            sid = h.split("|", 1)[0]
            if sid in seen:
                raise SystemExit(f"Duplicate sequence_id in marker FASTA: {sid}")
            seen.add(sid)
            s = clean_seq(s)
            if not s:
                continue
            out.write(f">{sid}\n{s}\n")
            n += 1
    print(f"prepared_sequences={n}")

if __name__ == "__main__":
    main()
