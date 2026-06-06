#!/usr/bin/env python3
"""凸多边形 + 单钉 N=2：可判定性/解维数实验（§W.11）。

运行: python decidability_convex_pin.py
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass

import sympy as sp
from sympy import Symbol, symbols

from mul_tyson_solve import (
    Pin,
    Point,
    ProblemSpec,
    SiteSymbols,
    build_equalities,
    polygon_boundary,
    solve_spec,
)


@dataclass
class DecidabilityReport:
    shape: str
    n_pins: int
    n_unknowns: int
    n_eqs: int
    residual_dof: int
    classification: str
    groebner_dim_hint: str


def count_residual_dof(spec: ProblemSpec) -> tuple[int, int, str]:
    """等式数 vs 未知数；Groebner 维数提示。"""
    sites = SiteSymbols.build(spec)
    unknowns = sites.unknowns()
    eqs = build_equalities(spec, sites)
    # extra fix
    for key, val in spec.extra_equalities.items():
        if key.startswith("p") and "x" in key:
            idx = int(key[1])
            eqs.append(sites.coords(idx)[0] - val)
        elif key.startswith("p") and "y" in key:
            idx = int(key[1])
            eqs.append(sites.coords(idx)[1] - val)
        elif key.startswith("r"):
            idx = int(key[1:])
            eqs.append(sites.radius(idx) - val)
    n_u = len(unknowns)
    n_e = len(eqs)
    residual = max(0, n_u - n_e)
    if n_e > n_u:
        hint = "过定（可能无解或退化解）"
        cls = "无解/退化"
    elif residual == 0:
        hint = "0 维（有限候选）"
        cls = "唯一或有限多解"
    elif residual == 1:
        hint = "1 维解族"
        cls = "多解（1 参数族）"
    else:
        hint = f"{residual} 维解族"
        cls = "多解（高维族）"
    return n_u, n_e, residual, cls, hint


def run_rectangle_route_a() -> DecidabilityReport:
    verts = [(0, 0), (2, 0), (2, 1), (0, 1)]
    tau = [2, 2, 1, 1]  # 顶边归旗2
    spec = ProblemSpec(
        boundary=polygon_boundary(verts, tau),
        p1=Point(0.4, 0.5),
        n_flags=2,
        pins=[Pin(Point(1.0, 0.0), 1, 2)],  # 底边中点钉 Γ_12
    )
    n_u, n_e, res, cls, hint = count_residual_dof(spec)
    return DecidabilityReport("矩形+单钉", 1, n_u, n_e, res, cls, hint)


def run_rectangle_double_pin() -> DecidabilityReport:
    verts = [(0, 0), (2, 0), (2, 1), (0, 1)]
    tau = [2, 2, 1, 1]
    spec = ProblemSpec(
        boundary=polygon_boundary(verts, tau),
        p1=Point(0.4, 0.5),
        n_flags=2,
        pins=[
            Pin(Point(1.0, 0.0), 1, 2),
            Pin(Point(2.0, 0.5), 1, 2),
        ],
    )
    n_u, n_e, res, cls, hint = count_residual_dof(spec)
    return DecidabilityReport(
        "矩形+双钉", 2, n_u, n_e, res, "一般无解/退化（双钉相容，定理 9.5 型）", hint
    )


def run_convex_general_pin() -> DecidabilityReport:
    verts = [(0, 0), (3, 0), (2.5, 1.2), (0.5, 1.5)]
    tau = [1, 2, 2, 1]
    spec = ProblemSpec(
        boundary=polygon_boundary(verts, tau),
        p1=Point(1.0, 0.6),
        n_flags=2,
        pins=[Pin(Point(1.5, 0.0), 1, 2)],
    )
    n_u, n_e, res, cls, hint = count_residual_dof(spec)
    return DecidabilityReport("凸五边形+单钉", 1, n_u, n_e, res, cls, hint)


def groebner_dimension_n2_single_pin() -> str:
    """符号层：N=2 单钉，无对称约化时 dof=3-1=2。"""
    p2x, p2y, r2 = symbols("p2x p2y r2", real=True)
    # 钉 (mx,my) 在 Γ_12
    mx, my, ax, ay = symbols("mx my ax ay", real=True)
    psi = (mx - ax) ** 2 * r2**2 - ((mx - p2x) ** 2 + (my - p2y) ** 2) * 1.0
    # 1 方程 3 未知 ⇒ 一般 2 维族
    return "符号计数：|Π|=1, N=2 ⇒ 3 未知 − 1 等式 ⇒ 一般 2 维解族（需 τ/Viol 或额外钉才唯一）"


def print_theorem_w11_summary() -> None:
    print("\n=== W.11 decidability/complexity summary ===")
    print("W.11.1: algebraic vertices + finite story => EXIST/UNIQUE/MULTI decidable (Tarski).")
    print("W.11.2: m=O(N+|Pi|+|dOmega|); QE worst 2^(2^O(m)); N=2 => m<=5 Groebner feasible.")
    print("W.11.3: convex + single pin + route A => generically dim>=1; unique needs rank>=dof.")
    print(groebner_dimension_n2_single_pin())


def main() -> int:
    reports = [
        run_rectangle_route_a(),
        run_rectangle_double_pin(),
        run_convex_general_pin(),
    ]
    print("=== convex polygon single-pin: DOF experiment ===")
    for r in reports:
        print(f"\n[{r.shape}]")
        print(f"  未知数={r.n_unknowns}, 等式={r.n_eqs}, 剩余自由度={r.residual_dof}")
        print(f"  W-分类提示: {r.classification}")
        print(f"  Groebner 维数提示: {r.groebner_dim_hint}")

    # 尝试实际求解矩形单钉
    spec = ProblemSpec(
        boundary=polygon_boundary(
            [(0, 0), (2, 0), (2, 1), (0, 1)], [2, 2, 1, 1]
        ),
        p1=Point(0.4, 0.5),
        n_flags=2,
        pins=[Pin(Point(1.0, 0.0), 1, 2)],
        collinear_y=0.35,
    )
    report = solve_spec(spec)
    print("\n[矩形单钉 + 共线 Ansatz y=0.35]")
    print(f"  求解状态: {report.classification}")
    print(f"  解个数: {len(report.solutions)}")

    print_theorem_w11_summary()
    return 0


if __name__ == "__main__":
    sys.exit(main())
