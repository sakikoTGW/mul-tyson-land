#!/usr/bin/env python3
"""Scan README.md for GitHub-unfriendly math patterns."""
import re
import sys
from pathlib import Path

README = Path(__file__).resolve().parents[1] / "README.md"
text = README.read_text(encoding="utf-8")
lines = text.splitlines()
issues = []

BANNED = [
    r"\\operatorname",
    r"\\dfrac",
    r"\\tag\{",
    r"\\label\{",
    r"\\text\{",
    r"\\require",
    r"\\html",
]

for i, line in enumerate(lines, 1):
    for pat in BANNED:
        if re.search(pat, line):
            issues.append((i, f"banned macro {pat}", line.strip()[:120]))

    # \frac without braces: \frac H or \frac x
    if re.search(r"\\frac\s+[^${\\]", line):
        issues.append((i, "frac missing braces", line.strip()[:120]))

    # table row: | inside $...$ (not \vert/\mid)
    if line.strip().startswith("|"):
        for m in re.finditer(r"\$([^$]+)\$", line):
            inner = m.group(1)
            if re.search(r"(?<!\\)\|", inner):
                issues.append((i, "raw | inside $ in table cell", m.group(0)[:80]))

# unbalanced $
for i, line in enumerate(lines, 1):
    if line.count("$") % 2:
        issues.append((i, "unbalanced $", line.strip()[:120]))

if issues:
    print(f"FOUND {len(issues)} issues:")
    for row in issues:
        print(f"  L{row[0]}: {row[1]} :: {row[2]}")
    sys.exit(1)
print("OK: no known GitHub math issues")
sys.exit(0)
