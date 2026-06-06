#!/usr/bin/env python3
"""Find README lines where GitHub inline $ math is likely to break."""
import re
import sys
from pathlib import Path

README = Path(__file__).resolve().parents[1] / "README.md"
lines = README.read_text(encoding="utf-8").splitlines()
issues = []

for i, line in enumerate(lines, 1):
    s = line.strip()
    if s.startswith("|") or s in ("$$",):
        continue
    if s.startswith("```"):
        continue
    n_dollar = line.count("$")
    if n_dollar < 2:
        continue
    # multiple inline segments or bold+underscore risk
    segs = re.findall(r"(?<!\$)\$([^$\n]+)\$(?!\$)", line)
    risk = len(segs) >= 2
    risk = risk or ("**" in line and segs)
    risk = risk or any("_" in seg or r"\*" in seg for seg in segs)
    if risk:
        issues.append((i, len(segs), line[:100]))

for row in issues:
    print(f"L{row[0]} ({row[1]} segs): {row[2]}")
print(f"\n{len(issues)} risky lines")
sys.exit(0)
