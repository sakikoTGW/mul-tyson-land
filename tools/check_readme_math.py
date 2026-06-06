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

RISKY_INLINE = re.compile(r"(?<!\$`)\$([^$`]+)\$(?!`)")

for i, line in enumerate(lines, 1):
    for pat in BANNED:
        if re.search(pat, line):
            issues.append((i, f"banned macro {pat}", line.strip()[:120]))

    if re.search(r"\\frac\s+[^${\\]", line):
        issues.append((i, "frac missing braces", line.strip()[:120]))

    if re.search(r"\*\*\$[^$]+\$|\$\*\*", line):
        issues.append((i, "bold wraps math", line.strip()[:120]))

    if re.search(r"dof\$", line):
        issues.append((i, "broken dof$", line.strip()[:120]))

    if re.search(r"边界归属：\$", line):
        issues.append((i, "old fragment boundary line", line.strip()[:120]))

    if re.search(r"\$`\\subset|\$`\\in V_1`\$", line):
        issues.append((i, "broken subset/in fragment", line.strip()[:120]))

    if "$$" not in line and re.search(r"(?<![\\lrvt])\|\\Omega\|", line):
        issues.append((i, "raw |Omega| outside display", line.strip()[:120]))

    if re.search(r"\$M\\in\\Gamma_\{12\}\$ [^\s(（]", line) and "仍在" not in line and "钉点 $M$" not in line:
        issues.append((i, "fragment $M\\in\\Gamma_{12}$", line.strip()[:120]))

    if line.strip().startswith("|"):
        for m in re.finditer(r"\$([^$]+)\$", line):
            inner = m.group(1)
            if re.search(r"(?<!\\)\|", inner):
                issues.append((i, "raw | inside $ in table", m.group(0)[:80]))

    # unprotected risky inline on label lines
    if (
        re.match(r"^\*\*定理 [0-9]", line)
        and "的证明" not in line
        and re.search(r"　\$[^$`]{10,}", line)
    ):
        issues.append((i, "theorem dense inline on label line", line.strip()[:120]))

    for m in RISKY_INLINE.finditer(line):
        inner = m.group(1)
        if re.fullmatch(r"\\Sigma_\d+", inner):
            continue
        if "_" in inner or r"\*" in inner:
            if not line.strip().startswith("```"):
                issues.append((i, "unprotected risky inline", m.group(0)[:80]))
                break

for i, line in enumerate(lines, 1):
    if line.count("$") % 2:
        issues.append((i, "unbalanced $", line.strip()[:120]))

if issues:
    print(f"FOUND {len(issues)} issues:")
    for row in issues[:40]:
        print(f"  L{row[0]}: {row[1]} :: {row[2]}")
    if len(issues) > 40:
        print(f"  ... and {len(issues)-40} more")
    sys.exit(1)
print("OK: no known GitHub math issues")
sys.exit(0)
