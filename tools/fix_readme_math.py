#!/usr/bin/env python3
"""Fix README math for GitHub KaTeX + table pipe safety."""
import re
from pathlib import Path

README = Path(__file__).resolve().parents[1] / "README.md"
text = README.read_text(encoding="utf-8")

# frac braces
text = re.sub(r"\\frac\s+H2\b", r"\\frac{H}{2}", text)
text = re.sub(r"\\frac\s+12\b", r"\\frac{1}{2}", text)
text = re.sub(r"\\frac\s+25\b", r"\\frac{2}{5}", text)
text = re.sub(r"\\frac\s+85\b", r"\\frac{8}{5}", text)

# cardinality / absolute value: |...| -> \lvert ... \rvert (inside math only)
def fix_bars_in_math(m):
    s = m.group(0)
    # eta|Omega| style
    s = re.sub(r"([^\\])\|(\\Omega\\b)", r"\1\\lvert\2", s)
    s = re.sub(r"(\\Omega)\\|", r"\1\\rvert", s)
    # |X| patterns
    s = re.sub(r"\|(\\mathfrak\{Sol\}[^$]*?)\\|", r"\\lvert\1\\rvert", s)
    s = re.sub(r"\|(\\mathcal\{E\})\\|", r"\\lvert\1\\rvert", s)
    s = re.sub(r"\|(\\Pi\\b)\\|", r"\\lvert\1\\rvert", s)
    s = re.sub(r"\|(\\mathfrak\{Sol\}\\cap\\mathfrak D\^c)\\|", r"\\lvert\1\\rvert", s)
    s = re.sub(r"\|(\\mathfrak\{Sol\})\\|", r"\\lvert\1\\rvert", s)
    return s

# Apply to each $...$ segment
def fix_segment(m):
    inner = m.group(1)
    # generic |expr| -> \lvert expr \rvert when expr doesn't contain |
    while True:
        new = re.sub(
            r"(?<![\\lrvt])\|([^|]+)\|",
            r"\\lvert\1\\rvert",
            inner,
        )
        if new == inner:
            break
        inner = new
    return "$" + inner + "$"

text = re.sub(r"\$([^$\n]+)\$", fix_segment, text)

README.write_text(text, encoding="utf-8")
print("fixed README.md")
