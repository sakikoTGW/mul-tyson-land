#!/usr/bin/env python3
"""扇形 Σ9 示意：路线 A 一参数族、份额 η 唯一、双钉 v1 无解区。"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Circle, Wedge


def sector_vertices(R: float, alpha: float, n: int = 48) -> list[tuple[float, float]]:
    pts = [(0.0, 0.0)]
    for k in range(n + 1):
        t = alpha * k / n
        pts.append((R * math.cos(t), R * math.sin(t)))
    return pts


def point_in_sector(x: float, y: float, R: float, alpha: float) -> bool:
    if x < -1e-12 or y < -1e-12:
        return False
    r = math.hypot(x, y)
    if r > R + 1e-12:
        return False
    return math.atan2(y, x) <= alpha + 1e-12


def assign_flag(x: float, y: float, a: float, b: float, lam: float) -> int:
    d1 = math.hypot(x - a, y)
    d2 = math.hypot(x - b, y) / lam
    return 1 if d1 <= d2 + 1e-9 else 2


def sector_area_v2(R: float, a: float, b: float, alpha: float, res: int = 220) -> float:
    lam = (R - b) / (R - a)
    cnt = 0
    total = 0
    for i in range(res):
        for j in range(res):
            x = R * i / (res - 1) * math.cos(alpha)
            y = R * j / (res - 1) * math.sin(alpha)
            if not point_in_sector(x, y, R, alpha):
                continue
            total += 1
            if assign_flag(x, y, a, b, lam) == 2:
                cnt += 1
    return cnt / max(total, 1)


def find_b_for_eta(R: float, a: float, alpha: float, eta: float) -> float:
    lo, hi = a + 1e-4, R - 1e-4
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if sector_area_v2(R, a, mid, alpha) > eta:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def plot_sector_case(
    R: float,
    a: float,
    alpha: float,
    b: float,
    title: str,
    out: Path,
    annotate_eta: float | None = None,
) -> None:
    lam = (R - b) / (R - a)
    p1, p2 = (a, 0.0), (b, 0.0)
    r1, r2 = 1.0, lam

    xs = np.linspace(-0.05, R + 0.15, 360)
    ys = np.linspace(-0.05, R * math.sin(alpha) + 0.15, 360)
    X, Y = np.meshgrid(xs, ys)
    labels = np.zeros_like(X)
    mask = np.zeros_like(X, dtype=bool)
    for iy in range(X.shape[0]):
        for ix in range(X.shape[1]):
            x, y = float(X[iy, ix]), float(Y[iy, ix])
            if point_in_sector(x, y, R, alpha):
                mask[iy, ix] = True
                labels[iy, ix] = assign_flag(x, y, a, b, lam)

    fig, ax = plt.subplots(figsize=(8.5, 7))
    c1, c2 = "#4C78A8", "#F58518"
    for flag, color in [(1, c1), (2, c2)]:
        region = np.ma.masked_where((labels != flag) | ~mask, labels)
        ax.contourf(X, Y, region, levels=[flag - 0.5, flag + 0.5], colors=[color], alpha=0.45)

    wedge = Wedge((0, 0), R, 0, math.degrees(alpha), fill=False, edgecolor="k", lw=2.2)
    ax.add_patch(wedge)
    ax.plot([0, R], [0, 0], "k-", lw=1.5)
    ax.plot([0, R * math.cos(alpha)], [0, R * math.sin(alpha)], "k-", lw=1.5)

    for cx, cy, r, c, name in [
        (p1[0], p1[1], r1, c1, "$p_1$"),
        (p2[0], p2[1], r2, c2, "$p_2$"),
    ]:
        ax.add_patch(Circle((cx, cy), r, fill=False, linestyle="--", edgecolor=c, lw=1.2))
        ax.scatter([cx], [cy], s=100, c=c, edgecolors="k", zorder=5)
        ax.annotate(name, (cx, cy), xytext=(6, 8), textcoords="offset points", fontsize=12)

    if annotate_eta is not None:
        share = sector_area_v2(R, a, b, alpha)
        ax.text(
            0.04,
            R * math.sin(alpha) * 0.88,
            rf"$\eta={annotate_eta:.2f}$  $\Rightarrow$  $b={b:.3f}$",
            fontsize=11,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
        )
        ax.text(
            0.04,
            R * math.sin(alpha) * 0.72,
            rf"$A(V_2)/|\Omega| \approx {share:.3f}$",
            fontsize=10,
        )

    ax.set_aspect("equal")
    ax.set_xlim(-0.05, R + 0.2)
    ax.set_ylim(-0.05, R * math.sin(alpha) + 0.2)
    ax.grid(True, alpha=0.25)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_v1_impossible(R: float, a: float, alpha: float, out: Path) -> None:
    """双钉 A、M：b 只能取 a 或 R²/a，均在 (a,R) 外或退化。"""
    b_roots = [a, R * R / a]
    fig, ax = plt.subplots(figsize=(8.5, 7))
    wedge = Wedge((0, 0), R, 0, math.degrees(alpha), fill=True, facecolor="#f0f0f0", edgecolor="k", lw=2.2)
    ax.add_patch(wedge)
    ax.scatter([a], [0], s=120, c="red", marker="x", linewidths=2, zorder=6, label="pin A")
    Mx = R * math.cos(alpha)
    ax.scatter([Mx], [0], s=120, c="darkred", marker="x", linewidths=2, zorder=6, label="pin M")
    for b in b_roots:
        color = "#888" if b > R or b < a else "#C44"
        ax.axvline(b, color=color, ls="--", lw=1.5)
        ax.text(b, R * math.sin(alpha) * 0.15, rf"$b={b:.2f}$", ha="center", fontsize=10)
    ax.axvspan(a, R, alpha=0.12, color="#4C78A8", label=r"feasible $b\in(a,R)$")
    ax.set_aspect("equal")
    ax.set_xlim(-0.05, R + 0.25)
    ax.set_ylim(-0.05, R * math.sin(alpha) + 0.15)
    ax.set_title(r"$\Sigma$9-I / v1: pins force $b\in\{a,\,R^2/a\}$ only", fontsize=11)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    out_dir = Path("figures")
    out_dir.mkdir(exist_ok=True)
    R, a, alpha = 2.0, 0.4, math.pi / 4

    plot_sector_case(
        R, a, alpha, b=1.1,
        title=r"Sector route A: sample $b=1.1$",
        out=out_dir / "fig04_sector_route_a.png",
    )
    eta = 0.35
    b_eta = find_b_for_eta(R, a, alpha, eta)
    plot_sector_case(
        R, a, alpha, b=b_eta,
        title=r"$\Sigma$9-IV: share $\eta$ fixes $b_\eta$ uniquely",
        out=out_dir / "fig05_sector_eta.png",
        annotate_eta=eta,
    )
    plot_v1_impossible(R, a, alpha, out_dir / "fig06_sector_v1.png")
    print("已保存扇形图:", out_dir.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
