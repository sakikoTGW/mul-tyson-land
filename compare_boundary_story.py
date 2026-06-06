#!/usr/bin/env python3
"""MWVD 正问题 / 容量泰森 / 边界故事反问题 — 系统对比（§W.9.10）。

不依赖 CGAL 安装；用同一乘性泰森模型数值演示三类工具的输入/输出差异。

运行: python compare_boundary_story.py [--out compare_report.txt]
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from mul_tyson_solve import (
    Pin,
    Point,
    ProblemSpec,
    assign_min_flag,
    polygon_boundary,
    solve_spec,
    verify_viol,
)


# ---------------------------------------------------------------------------
# 正问题：给定 (p_i, r_i) 画 MWVD（与 CGAL Apollonius / ESRI 同型）
# ---------------------------------------------------------------------------


def mwvd_assign(
    x: float, y: float, sites: list[tuple[float, float]], radii: list[float]
) -> int:
    return assign_min_flag(x, y, sites, radii)


def edge_violation_rate(
    vertices: list[tuple[float, float]],
    tau: list[int],
    sites: list[tuple[float, float]],
    radii: list[float],
    samples: int = 40,
) -> tuple[float, list[tuple[int, int, float]]]:
    """边界故事违约率：σ(s)≠τ(s) 的边比例。"""
    n = len(vertices)
    bad_edges: list[tuple[int, int, float]] = []
    for k in range(n):
        x0, y0 = vertices[k]
        x1, y1 = vertices[(k + 1) % n]
        mism = 0
        for t in range(samples):
            u = t / (samples - 1) if samples > 1 else 0.0
            x, y = x0 + u * (x1 - x0), y0 + u * (y1 - y0)
            sigma = mwvd_assign(x, y, sites, radii)
            if sigma != tau[k]:
                mism += 1
        rate = mism / samples
        if rate > 0.05:
            bad_edges.append((k, tau[k], rate))
    total_bad = len(bad_edges)
    return total_bad / n, bad_edges


def cell_area_ratio(
    vertices: list[tuple[float, float]],
    sites: list[tuple[float, float]],
    radii: list[float],
    flag: int,
    grid: int = 120,
) -> float:
    xs = [v[0] for v in vertices]
    ys = [v[1] for v in vertices]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    cnt = 0
    total = 0
    for i in range(grid):
        for j in range(grid):
            x = xmin + (xmax - xmin) * (i + 0.5) / grid
            y = ymin + (ymax - ymin) * (j + 0.5) / grid
            # 粗判：在 bbox 内（凸矩形算例足够）
            if not (xmin <= x <= xmax and ymin <= y <= ymax):
                continue
            total += 1
            if mwvd_assign(x, y, sites, radii) == flag:
                cnt += 1
    return cnt / max(total, 1)


# ---------------------------------------------------------------------------
# 容量泰森：迭代调 r2 使 Area(V2)/Area(Ω) ≈ η（Balzer ISVD 2009 型）
# ---------------------------------------------------------------------------


def capacity_adjust_r2(
    vertices: list[tuple[float, float]],
    p1: tuple[float, float],
    p2: tuple[float, float],
    r1: float,
    target_eta: float,
    r2_init: float = 0.5,
    max_iter: int = 40,
) -> float:
    """只调权重 r2，不移动站点 — 标准容量泰森启发式。"""
    r2 = r2_init
    for _ in range(max_iter):
        ratio = cell_area_ratio(vertices, [p1, p2], [r1, r2], flag=2)
        if abs(ratio - target_eta) < 0.01:
            break
        # 面积偏大 ⇒ 增大 r2（减弱旗2）
        if ratio > target_eta:
            r2 *= 1.08
        else:
            r2 *= 0.92
        r2 = max(0.05, min(r2, 5.0))
    return r2


# ---------------------------------------------------------------------------
# 对比场景
# ---------------------------------------------------------------------------


@dataclass
class ScenarioResult:
    name: str
    method: str
    p2: tuple[float, float] | None
    r2: float | None
    viol_rate: float
    bad_edges: list
    area_ratio_v2: float
    notes: str = ""


@dataclass
class ComparisonReport:
    scenario: str
    tau: list[int]
    results: list[ScenarioResult] = field(default_factory=list)


def _solve_r2_from_pin(
    mx: float,
    my: float,
    p1: tuple[float, float],
    p2x: float,
    p2y: float,
) -> float | None:
    """由 Psi_12(M)=0 解 r2（正根）。"""
    ax, ay = p1
    d1 = (mx - ax) ** 2 + (my - ay) ** 2
    d2 = (mx - p2x) ** 2 + (my - p2y) ** 2
    if d2 <= 1e-12:
        return None
    r2 = math.sqrt(d1 / d2)
    return r2 if r2 > 1e-6 else None


def _best_inverse_solution(
    spec: ProblemSpec, sol, p1: tuple[float, float]
) -> tuple[tuple[float, float] | None, float | None, bool]:
    """从有限解或 1 维族扫描 Viol PASS 的解。"""
    if sol.solutions:
        s = sol.solutions[0]
        p2 = (float(s["p2x"]), float(s.get("p2y", spec.collinear_y or 0.0)))
        r2 = float(s["r2"])
        ok, _ = verify_viol(spec, [p1, p2], [1.0, r2])
        return p2, r2, ok
    if not spec.pins:
        return None, None, False
    pin = spec.pins[0]
    p2y = spec.collinear_y if spec.collinear_y is not None else 0.35
    for k in range(80):
        p2x = 0.55 + 1.35 * k / 79
        r2 = _solve_r2_from_pin(pin.point.x, pin.point.y, p1, p2x, p2y)
        if r2 is None:
            continue
        p2 = (p2x, p2y)
        ok, _ = verify_viol(spec, [p1, p2], [1.0, r2], samples_per_edge=16)
        if ok:
            return p2, r2, True
    return None, None, False


def run_rectangle_scenario() -> ComparisonReport:
    """矩形：顶边要归旗2，底边中点钉 Γ_12 — 边界故事驱动。"""
    verts = [(0, 0), (2, 0), (2, 1), (0, 1)]
    tau = [2, 2, 1, 1]  # 边0底、边1右、边2顶→2、边3左→1
    p1 = (0.4, 0.5)
    rep = ComparisonReport("矩形：顶边归旗2 + 底中点钉", tau)

    # 1) CGAL/ESRI 型：用户随便给一组参数
    p2_guess = (1.2, 0.35)
    r2_guess = 0.6
    vr, be = edge_violation_rate(verts, tau, [p1, p2_guess], [1.0, r2_guess])
    ar = cell_area_ratio(verts, [p1, p2_guess], [1.0, r2_guess], 2)
    rep.results.append(
        ScenarioResult(
            "MWVD正问题(手选参数)",
            "CGAL/ESRI 同型",
            p2_guess,
            r2_guess,
            vr,
            be,
            ar,
            "正问题：给定参数即画划分；不保证 τ",
        )
    )

    # 2) 容量泰森：调 r2 使面积≈30%，站点不动
    target_eta = 0.30
    r2_cap = capacity_adjust_r2(verts, p1, p2_guess, 1.0, target_eta)
    vr2, be2 = edge_violation_rate(verts, tau, [p1, p2_guess], [1.0, r2_cap])
    ar2 = cell_area_ratio(verts, [p1, p2_guess], [1.0, r2_cap], 2)
    rep.results.append(
        ScenarioResult(
            "容量泰森",
            "Balzer ISVD 2009 型",
            p2_guess,
            r2_cap,
            vr2,
            be2,
            ar2,
            f"面积份额≈{target_eta:.0%}；仍不保证顶边归旗2",
        )
    )

    # 3) 本文：反求 (p2,r2) 满足钉点 + Viol
    spec = ProblemSpec(
        boundary=polygon_boundary(verts, tau),
        p1=Point(*p1),
        n_flags=2,
        pins=[Pin(Point(1.0, 0.0), 1, 2)],
        collinear_y=0.35,
    )
    sol = solve_spec(spec)
    p2_inv, r2_inv, ok = _best_inverse_solution(spec, sol, p1)
    if p2_inv is not None and r2_inv is not None:
        vr3, be3 = edge_violation_rate(verts, tau, [p1, p2_inv], [1.0, r2_inv])
        ar3 = cell_area_ratio(verts, [p1, p2_inv], [1.0, r2_inv], 2)
        rep.results.append(
            ScenarioResult(
                "边界故事反问题",
                "本文 Viol+Π",
                p2_inv,
                r2_inv,
                vr3,
                be3,
                ar3,
                f"Viol={'PASS' if ok else 'FAIL'}; class={sol.classification}",
            )
        )
    else:
        rep.results.append(
            ScenarioResult(
                "边界故事反问题",
                "本文 Viol+Π",
                None,
                None,
                1.0,
                [],
                0.0,
                f"求解失败: {sol.classification}",
            )
        )

    return rep


def format_report(rep: ComparisonReport) -> str:
    lines = [
        "=" * 72,
        f"场景: {rep.scenario}",
        f"边界故事 τ (按边): {rep.tau}",
        "=" * 72,
        "",
        "| 方法 | p2 | r2 | 违约边比例 | V2面积比 | 备注 |",
        "|------|-----|-----|------------|----------|------|",
    ]
    for r in rep.results:
        p2s = f"({r.p2[0]:.3f},{r.p2[1]:.3f})" if r.p2 else "—"
        r2s = f"{r.r2:.4f}" if r.r2 is not None else "—"
        lines.append(
            f"| {r.method} | {p2s} | {r2s} | {r.viol_rate:.0%} ({len(r.bad_edges)}边) | "
            f"{r.area_ratio_v2:.0%} | {r.notes} |"
        )
    lines.extend(
        [
            "",
            "结论:",
            "  - MWVD forward (CGAL/ESRI): sites in, partition out; no boundary story tau.",
            "  - Capacity Voronoi: tune weights for area; moment constraints, not edge labels.",
            "  - This work: inverse solve for Viol=0 + pins Pi; three different problems.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    reports = [run_rectangle_scenario()]
    text = "\n".join(format_report(r) for r in reports)
    print(text)
    if args.out:
        args.out.write_text(text, encoding="utf-8")
        print(f"已写入 {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
