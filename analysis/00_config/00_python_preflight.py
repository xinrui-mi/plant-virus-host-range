#!/usr/bin/env python3
import sys, importlib.util
mods = ["Bio", "numpy"]
missing = [m for m in mods if importlib.util.find_spec(m) is None]
if missing:
    print("Missing Python modules:", ", ".join(missing), file=sys.stderr)
    sys.exit(2)
print("Python preflight PASS: Bio, numpy")
