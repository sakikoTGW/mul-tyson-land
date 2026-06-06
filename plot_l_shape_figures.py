#!/usr/bin/env python3
"""L 形示意：§15.1 两面旗月牙 + §21.2 四旗路线 C（论文标准参数）。"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

from mul_tyson_solve import Pin, Point, ProblemSpec, polygon_boundary
from mul_tyson_viz import DomainInput, InteriorCell, plot_tyson

# §15.0 标准 L：R_x=2, w=0.5, R_y=2
L_VERTICES = [(0, 0), (2, 0), (2, 0.5), (0.5, 0.5), (0.5, 2), (0, 2)]
RX, W, RY = 2.0, 0.5, 2.0
A1, B1 = 0.4, 0.6
MID = W / 2  # 0.25
LAM = (RX - B1) / (RX - A1)  # 7/8


def l_domain() -> DomainInput:
    return DomainInput(vertices=L_VERTICES, holes=[])


def _spec_two_flag() -> ProblemSpec:
    tau = [1, 2, 1, 1, 1, 1]
    pin = Pin(Point(RX, MID), 1, 2)  # M_x
    p1 = Point(A1, MID)
    return ProblemSpec(
        boundary=polygon_boundary(L_VERTICES, tau),
        p1=p1,
        n_flags=2,
        pins=[pin],
        interior=[
            InteriorCell(Point(0.45, MID), 1),
            InteriorCell(Point(1.2, MID), 2),
        ],
        collinear_y=MID,
    )


def _spec_four_flag() -> ProblemSpec:
    xi, eps, r = 11 / 20, 1 / 10, 9 / 10
    tau = [1, 2, 1, 1, 1, 1]
    pin = Pin(Point(RX, MID), 1, 2)
    p1 = Point(A1, MID)
    return ProblemSpec(
        boundary=polygon_boundary(L_VERTICES, tau),
        p1=p1,
        n_flags=4,
        pins=[pin],
        interior=[
            InteriorCell(Point(0.45, MID), 1),
            InteriorCell(Point(1.2, MID), 2),
            InteriorCell(Point(xi, W - eps - 0.02), 3),
            InteriorCell(Point(xi, eps + 0.02), 4),
        ],
        collinear_y=MID,
    )


def plot_two_flag(out: Path) -> None:
    spec = _spec_two_flag()
    sol = {"p2x": B1, "p2y": MID, "r2": LAM}
    plot_tyson(l_domain(), spec, sol, out, dpi=150, grid_res=400)
    print(f"§15.1 两面旗: p2=({B1},{MID}), λ={LAM:.4f} -> {out}")


def plot_four_flag(out: Path) -> None:
    xi, eps, r = 11 / 20, 1 / 10, 9 / 10
    spec = _spec_four_flag()
    sol = {
        "p2x": B1,
        "p2y": MID,
        "r2": LAM,
        "p3x": xi,
        "p3y": W - eps,
        "r3": r,
        "p4x": xi,
        "p4y": eps,
        "r4": r,
    }
    plot_tyson(l_domain(), spec, sol, out, dpi=150, grid_res=400)
    print(f"§21.2 四旗路线C: p3/p4=({xi:.2f},{W-eps:.2f})/({xi:.2f},{eps:.2f}), r={r} -> {out}")


def plot_crescent_share(out: Path) -> None:
    """水平臂内 V2 面积占比（命题 15.3 ~73%）。"""
    p1, p2 = (A1, MID), (B1, MID)
    res = 300
    cnt2 = cnt = 0
    for i in range(res):
        for j in range(res):
            x = RX * i / (res - 1)
            y = W * j / (res - 1)
            cnt += 1
            d1 = math.hypot(x - p1[0], y - p1[1])
            d2 = math.hypot(x - p2[0], y - p2[1]) / LAM
            if d2 < d1 - 1e-9:
                cnt2 += 1
    share = cnt2 / max(cnt, 1)

    spec = _spec_two_flag()
    sol = {"p2x": B1, "p2y": MID, "r2": LAM}
    from mul_tyson_viz import compute_region_grid, bbox

    pos = [(A1, MID), (B1, MID)]
    rad = [1.0, LAM]
    X, Y, labels, mask = compute_region_grid(L_VERTICES, [], pos, rad, res=400)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for flag, color in [(1, "#4C78A8"), (2, "#F58518")]:
        region = np.ma.masked_where((labels != flag) | ~mask, labels)
        ax.contourf(X, Y, region, levels=[flag - 0.5, flag + 0.5], colors=[color], alpha=0.5)
    ox, oy = zip(*L_VERTICES)
    ax.plot(list(ox) + [ox[0]], list(oy) + [oy[0]], "k-", lw=2)
    ax.axvline(W, color="k", ls=":", lw=1, alpha=0.5)
    for cx, cy, r, c, name in [
        (A1, MID, 1.0, "#4C78A8", "$p_1$"),
        (B1, MID, LAM, "#F58518", "$p_2$"),
    ]:
        ax.add_patch(Circle((cx, cy), r, fill=False, ls="--", edgecolor=c, lw=1.1))
        ax.scatter([cx], [cy], s=90, c=c, edgecolors="k", zorder=5)
        ax.annotate(name, (cx, cy), xytext=(6, 6), textcoords="offset points", fontsize=11)
    ax.scatter([RX], [MID], marker="x", c="red", s=70, zorder=6)
    ax.annotate("$M_x$", (RX, MID), xytext=(-28, 8), textcoords="offset points", color="red")
    ax.set_aspect("equal")
    xmin, xmax, ymin, ymax = bbox(L_VERTICES, 0.05)
    ax.set_xlim(xmin, xmax + 0.1)
    ax.set_ylim(ymin, ymax + 0.1)
    ax.set_title(
        rf"L horizontal arm: $A(V_2)/|\mathrm{{arm}}| \approx {share:.3f}$",
        fontsize=11,
    )
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"月牙占比 ~{share:.1%} -> {out}")


def main() -> int:
    out_dir = Path("figures")
    out_dir.mkdir(exist_ok=True)
    plot_two_flag(out_dir / "fig10_L_tau.png")
    plot_two_flag(out_dir / "fig02_L_auto.png")
    plot_four_flag(out_dir / "fig11_L_routeC.png")
    plot_crescent_share(out_dir / "fig12_L_crescent.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
