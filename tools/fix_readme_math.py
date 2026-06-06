#!/usr/bin/env python3
"""Fix README math for GitHub KaTeX rendering."""
import re
from pathlib import Path

README = Path(__file__).resolve().parents[1] / "README.md"
text = README.read_text(encoding="utf-8")

# --- structural prose fixes (idempotent) ---
REPL = [
    ("$M\\in\\Gamma_{12}$ 不变", "钉点 $M$ 仍在 $\\Gamma_{12}$ 上"),
    ("钉点 $M\\in\\Gamma_{12}$ 仍满足", "钉点 $M$ 仍满足"),
    ("而不移动 $M\\in\\Gamma_{12}$", "而钉点 $M$ 不移动"),
    ("未知 $(b,\\lambda)$；钉 $M\\in\\Gamma_{12}$ 得", "未知 $(b,\\lambda)$；钉点 $M$ 在 $\\Gamma_{12}$ 上，得"),
    ("**定理 $\\Sigma 9$**", "定理 $\\Sigma 9$"),
    (
        "边界归属：$\\sigma_{\\mathcal S}(s)=\\tau(s)$ 对 $s\\in\\partial\\Omega$ 成立（分界点除外）。其中",
        "在边界 $\\partial\\Omega$ 上（分界点除外），还要求归属一致：",
    ),
    ("dof$=3$", r"$\mathrm{dof}=3$"),
    (r"\mathcal B^\*", r"\mathcal B^{\ast}"),
    (r"\arg\min_i", r"\arg\min_{i}"),
    (r"\left\{", r"\\{"),
    (r"\right\}", r"\\}"),
]
for old, new in REPL:
    text = text.replace(old, new)

text = re.sub(r"\\frac\s+H2\b", r"\\frac{H}{2}", text)
text = re.sub(r"\*\*\$([^$]+)\$\*\*", r"$\1$", text)
text = re.sub(r"\*\*(\$\\Sigma 9\$)-([IV]+)\*\*", r"\1-\2", text)

# --- intro: display blocks for theorems (GitHub breaks dense inline on ** lines) ---
INTRO_OLD = """设 $\\Omega\\subset\\mathbb R^2$ 为有界闭域，$\\partial\\Omega$ 由分段 $C^1$ 折线围成；$\\Omega$ 的边界可为任意多边形（含海星形等非凸外形），亦允许孔洞 $\\partial\\Omega=\\Gamma_{\\mathrm{out}}\\sqcup\\bigsqcup_h\\Gamma_{\\mathrm{in}}^{(h)}$（例如环域：外圈为 $\\Gamma_{\\mathrm{out}}$，内圈为孔边界 $\\Gamma_{\\mathrm{in}}^{(1)}$）。给定 $p_1\\in\\Omega$，第一面旗半径归一 $r_1=1$。对其余旗 $(p_i,r_i)$（$r_i>0$），定义

$$
\\varphi_i(x):=\\frac{\\|x-p_i\\|}{r_i},
$$

势力区 $V_i=\\{x\\in\\Omega:\\varphi_i(x)\\le\\varphi_j(x),\\,\\forall j\\neq i\\}$。旗 $i$ 与 $j$ 的分界线 $\\{\\varphi_i=\\varphi_j\\}$ 为阿波罗尼奥斯圆；$r_i=r_j$ 时退化为直线。边界标签 $\\tau$ 指定 $\\partial\\Omega$ 各边（除分界点外）的目标归属。要求

$$
\\bigcup_{i=1}^N V_i=\\Omega.
$$

在边界 $\\partial\\Omega$ 上（分界点除外），还要求归属一致：

$$
\\sigma_{\\mathcal S}(s)=\\tau(s),\\qquad s\\in\\partial\\Omega,
$$

其中

$$
\\sigma_{\\mathcal S}(x)=\\arg\\min_{i}\\varphi_i(x).
$$

**反问题**：已知 $\\Omega,p_1,\\tau$，求 $(p_2,\\ldots,p_N,r_2,\\ldots,r_N)$，并判定 $\\mathfrak{Sol}$ 为无解、唯一、有限或无穷；无穷时给出参数化闭式。

在路线 A（(1.2)(1.3)）与强规范 $\\mathcal B^{\\ast}$（渗漏边集 $\\mathcal E=\\varnothing$，见 §4）下，记 $m_{\\mathrm{eff}}=\\lvert\\mathcal E\\rvert$，有

**定理 1.1**　$m_{\\mathrm{eff}}=0\\Rightarrow n_{\\min}=2$；$m_{\\mathrm{eff}}\\ge1\\Rightarrow n_{\\min}=3$，且可取 $(p_3,r_3)=(p_2,\\lambda)$。

**定理 1.2**　$\\Omega=[0,W]\\times[0,H]$，$p_1=(a,H/2)$，$M=(W,H/2)\\in\\Gamma_{12}$ 时，$b=W-(W-a)r_2$（$r_2>0$）为 $\\mathfrak{Sol}$ 的一条 1 维族。

**定理 1.3**　扇形 $\\Omega(R,\\alpha)$、路线 A 的 1 维族上，若 $\\mathrm{Area}(V_2)=\\eta\\lvert\\Omega\\rvert$（$\\eta\\in(0,1)$），则存在唯一 $b_\\eta\\in(a,R)$。"""

INTRO_NEW = """设 $\\Omega\\subset\\mathbb R^2$ 为有界闭域，$\\partial\\Omega$ 由分段 $C^1$ 折线围成。$\\Omega$ 的边界可为任意多边形（含海星形等非凸外形）；亦允许孔洞，例如环域

$$
\\partial\\Omega=\\Gamma_{\\mathrm{out}}\\sqcup\\bigsqcup_h\\Gamma_{\\mathrm{in}}^{(h)}.
$$

外圈为 $\\Gamma_{\\mathrm{out}}$，内圈为孔边界 $\\Gamma_{\\mathrm{in}}^{(h)}$。给定 $p_1\\in\\Omega$，第一面旗半径归一 $r_1=1$。对其余旗 $(p_i,r_i)$（$r_i>0$），定义

$$
\\varphi_i(x):=\\frac{\\|x-p_i\\|}{r_i},
$$

势力区定义为

$$
V_i=\\{x\\in\\Omega:\\varphi_i(x)\\le\\varphi_j(x),\\,\\forall j\\neq i\\}.
$$

旗 $i$ 与 $j$ 的分界线 $\\{\\varphi_i=\\varphi_j\\}$ 为阿波罗尼奥斯圆；$r_i=r_j$ 时退化为直线。边界标签 $\\tau$ 指定 $\\partial\\Omega$ 各边（除分界点外）的目标归属。要求

$$
\\bigcup_{i=1}^N V_i=\\Omega.
$$

在边界 $\\partial\\Omega$ 上（分界点除外），还要求归属一致：

$$
\\sigma_{\\mathcal S}(s)=\\tau(s),\\qquad s\\in\\partial\\Omega,
$$

其中

$$
\\sigma_{\\mathcal S}(x)=\\arg\\min_{i}\\varphi_i(x).
$$

**反问题**　已知 $\\Omega,p_1,\\tau$，求其余旗位与半径，并判定解集为无解、唯一、有限或无穷；无穷时给出参数化闭式。

在路线 A（(1.2)(1.3)）与强规范 $\\mathcal B^{\\ast}$（渗漏边集 $\\mathcal E=\\varnothing$，见 §4）下，记

$$
m_{\\mathrm{eff}}=\\lvert\\mathcal E\\rvert.
$$

**定理 1.1**

$$
m_{\\mathrm{eff}}=0\\Rightarrow n_{\\min}=2;\\qquad m_{\\mathrm{eff}}\\ge1\\Rightarrow n_{\\min}=3,\\quad (p_3,r_3)=(p_2,\\lambda).
$$

**定理 1.2**

$$
\\Omega=[0,W]\\times[0,H],\\quad p_1=(a,H/2),\\quad M=(W,H/2)\\in\\Gamma_{12}
\\ \\Longrightarrow\\ b=W-(W-a)r_2\\quad (r_2>0),
$$

为解集的一条 1 维族。

**定理 1.3**　扇形 $\\Omega(R,\\alpha)$、路线 A 的 1 维族上，若

$$
\\mathrm{Area}(V_2)=\\eta\\lvert\\Omega\\rvert,\\qquad \\eta\\in(0,1),
$$

则存在唯一 $b_\\eta\\in(a,R)$。"""

# apply intro if still old form
if "**定理 1.1**　$m" in text or "边界归属：" in text:
    if "边界归属：" in text:
        pass  # partial; do piecemeal below
    if "**定理 1.1**　$m" in text:
        start = text.find("设 $\\Omega")
        end = text.find("---\n\n## 2")
        if start != -1 and end != -1:
            # rebuild intro from INTRO_NEW - match from file start markers
            pass

# piecemeal intro if full block replace fails
blocks = [
    (
        "设 $\\Omega\\subset\\mathbb R^2$ 为有界闭域，$\\partial\\Omega$ 由分段 $C^1$ 折线围成；$\\Omega$ 的边界可为任意多边形（含海星形等非凸外形），亦允许孔洞 $\\partial\\Omega=\\Gamma_{\\mathrm{out}}\\sqcup\\bigsqcup_h\\Gamma_{\\mathrm{in}}^{(h)}$（例如环域：外圈为 $\\Gamma_{\\mathrm{out}}$，内圈为孔边界 $\\Gamma_{\\mathrm{in}}^{(1)}$）。给定 $p_1\\in\\Omega$，第一面旗半径归一 $r_1=1$。对其余旗 $(p_i,r_i)$（$r_i>0$），定义",
        "设 $\\Omega\\subset\\mathbb R^2$ 为有界闭域，$\\partial\\Omega$ 由分段 $C^1$ 折线围成。$\\Omega$ 的边界可为任意多边形（含海星形等非凸外形）；亦允许孔洞，例如环域\n\n$$\n\\partial\\Omega=\\Gamma_{\\mathrm{out}}\\sqcup\\bigsqcup_h\\Gamma_{\\mathrm{in}}^{(h)}.\n$$\n\n外圈为 $\\Gamma_{\\mathrm{out}}$，内圈为孔边界 $\\Gamma_{\\mathrm{in}}^{(h)}$。给定 $p_1\\in\\Omega$，第一面旗半径归一 $r_1=1$。对其余旗 $(p_i,r_i)$（$r_i>0$），定义",
    ),
    (
        "势力区 $V_i=\\{x\\in\\Omega:\\varphi_i(x)\\le\\varphi_j(x),\\,\\forall j\\neq i\\}$。旗 $i$ 与 $j$ 的分界线 $\\{\\varphi_i=\\varphi_j\\}$ 为阿波罗尼奥斯圆；$r_i=r_j$ 时退化为直线。边界标签 $\\tau$ 指定 $\\partial\\Omega$ 各边（除分界点外）的目标归属。要求",
        "势力区定义为\n\n$$\nV_i=\\{x\\in\\Omega:\\varphi_i(x)\\le\\varphi_j(x),\\,\\forall j\\neq i\\}.\n$$\n\n旗 $i$ 与 $j$ 的分界线 $\\{\\varphi_i=\\varphi_j\\}$ 为阿波罗尼奥斯圆；$r_i=r_j$ 时退化为直线。边界标签 $\\tau$ 指定 $\\partial\\Omega$ 各边（除分界点外）的目标归属。要求",
    ),
    (
        "**定理 1.1**　$m_{\\mathrm{eff}}=0\\Rightarrow n_{\\min}=2$；$m_{\\mathrm{eff}}\\ge1\\Rightarrow n_{\\min}=3$，且可取 $(p_3,r_3)=(p_2,\\lambda)$。",
        "**定理 1.1**\n\n$$\nm_{\\mathrm{eff}}=0\\Rightarrow n_{\\min}=2;\\qquad m_{\\mathrm{eff}}\\ge1\\Rightarrow n_{\\min}=3,\\quad (p_3,r_3)=(p_2,\\lambda).\n$$",
    ),
    (
        "**定理 1.2**　$\\Omega=[0,W]\\times[0,H]$，$p_1=(a,H/2)$，$M=(W,H/2)\\in\\Gamma_{12}$ 时，$b=W-(W-a)r_2$（$r_2>0$）为 $\\mathfrak{Sol}$ 的一条 1 维族。",
        "**定理 1.2**\n\n$$\n\\Omega=[0,W]\\times[0,H],\\quad p_1=(a,H/2),\\quad M=(W,H/2)\\in\\Gamma_{12}\n\\ \\Longrightarrow\\ b=W-(W-a)r_2\\quad (r_2>0),\n$$\n\n为解集的一条 1 维族。",
    ),
    (
        "**定理 1.3**　扇形 $\\Omega(R,\\alpha)$、路线 A 的 1 维族上，若 $\\mathrm{Area}(V_2)=\\eta\\lvert\\Omega\\rvert$（$\\eta\\in(0,1)$），则存在唯一 $b_\\eta\\in(a,R)$。",
        "**定理 1.3**　扇形 $\\Omega(R,\\alpha)$、路线 A 的 1 维族上，若\n\n$$\n\\mathrm{Area}(V_2)=\\eta\\lvert\\Omega\\rvert,\\qquad \\eta\\in(0,1),\n$$\n\n则存在唯一 $b_\\eta\\in(a,R)$。",
    ),
    (
        "在路线 A（(1.2)(1.3)）与强规范 $\\mathcal B^{\\ast}$（渗漏边集 $\\mathcal E=\\varnothing$，见 §4）下，记 $m_{\\mathrm{eff}}=\\lvert\\mathcal E\\rvert$，有",
        "在路线 A（(1.2)(1.3)）与强规范 $\\mathcal B^{\\ast}$（渗漏边集 $\\mathcal E=\\varnothing$，见 §4）下，记\n\n$$\nm_{\\mathrm{eff}}=\\lvert\\mathcal E\\rvert.\n$$",
    ),
    (
        "- 允许**孔洞**：$\\partial\\Omega=\\Gamma_{\\mathrm{out}}\\sqcup\\bigsqcup_h\\Gamma_{\\mathrm{in}}^{(h)}$，孔边界与外边界的归属均由 $\\tau$ 指定；",
        "- 允许**孔洞**。边界分解为\n\n  $$\n  \\partial\\Omega=\\Gamma_{\\mathrm{out}}\\sqcup\\bigsqcup_h\\Gamma_{\\mathrm{in}}^{(h)}.\n  $$\n\n  孔边界与外边界的归属均由 $\\tau$ 指定；",
    ),
]
for old, new in blocks:
    text = text.replace(old, new)

# GitHub: wrap risky inline math with backticks inside $...$
RISKY = re.compile(r"_|\\\*|\\mathfrak|\\mathrm|\\mathcal|\\varnothing|\\Rightarrow|\\Longleftrightarrow|\\boxed|\\lvert|\\rvert|\\sqcup|\\bigsqcup|\\varphi|\\Gamma|\\Omega|\\partial|\\mathcal|\\mathfrak|\\dim|\\min|\\max|\\eta|\\lambda|\\chi|\\tau|\\sigma|\\varnothing")

def protect_segment(m):
    inner = m.group(1)
    if inner.startswith("`") and inner.endswith("`"):
        return m.group(0)
    if RISKY.search(inner):
        return "$`" + inner + "`$"
    return "$" + inner + "$"

def protect_non_display(chunk):
    return re.sub(r"\$([^$\n]+)\$", protect_segment, chunk)

def fix_display_block(part: str) -> str:
    inner = part[2:-2] if part.startswith("$$") and part.endswith("$$") else part
    inner = re.sub(r"(?<![\\lrvt])\|(\\Omega)\|", r"\\lvert\1\\rvert", inner)
    inner = re.sub(r"(?<![\\lrvt])\|(\\varphi)\|", r"\\lvert\1\\rvert", inner)
    return "$$" + inner + "$$" if part.startswith("$$") else inner

FRAGMENT_REPL = [
    ("$`\\subset V_1`$", r"$\subset V_1$"),
    ("邻域 $`\\in V_1`$", "邻域 $\\subset V_1$"),
    ("$`\\subset V_1`$", r"$\subset V_1$"),
]
for old, new in FRAGMENT_REPL:
    text = text.replace(old, new)

parts = re.split(r"(\$\$.*?\$\$)", text, flags=re.DOTALL)
for i, part in enumerate(parts):
    if part.startswith("$$"):
        parts[i] = fix_display_block(part)
    else:
        parts[i] = protect_non_display(part)
text = "".join(parts)

# table pipes inside math
def fix_segment(m):
    inner = m.group(1)
    while True:
        new = re.sub(r"(?<![\\lrvt])\|([^|]+)\|", r"\\lvert\1\\rvert", inner)
        if new == inner:
            break
        inner = new
    return "$" + inner + "$"

# only for non-backtick segments in tables - skip if already has backticks
def fix_tables(chunk):
    lines = chunk.split("\n")
    out = []
    for line in lines:
        if line.strip().startswith("|") and "$`" not in line:
            line = re.sub(r"\$([^$\n]+)\$", fix_segment, line)
        out.append(line)
    return "\n".join(out)

text = fix_tables(text)

README.write_text(text, encoding="utf-8")
print("fixed README.md")
