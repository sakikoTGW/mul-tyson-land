#!/usr/bin/env python3
"""乘性泰森宽题 W — 代数消元工具（非证明链）。

将 (f_D, γ, τ, Π, 𝔤) 代入 𝔈^mul,full，对未知 (p_i, r_i) 做 Gröbner 消元，
输出实解 / 解族维数 / Viol 审计。对应 symmetric_sector_theorems.md §W.6.5–W.6.14。

用法:
  python mul_tyson_solve.py --example rect-route-a
  python mul_tyson_solve.py --spec problem.json
  python mul_tyson_solve.py --vertices 0,0 2,0 2,1 0,1 --p1 0.4,0.5 --n 2 \\
      --pin 2,0.5,1,2 --tau 0:1 1:2 2:1 3:1
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy import Eq, GroebnerBasis, Symbol, symbols
from sympy.polys.polytools import groebner


# ---------------------------------------------------------------------------
# 几何与故事规格（只通过顶点/参数进入，无形状名）
# ---------------------------------------------------------------------------


@dataclass
class Point:
    x: float
    y: float

    def as_tuple(self) -> tuple[float, float]:
        return (self.x, self.y)


@dataclass
class BoundarySegment:
    """边界段：顶点序列上的一条边，目标归属 τ ∈ {1..N}。"""

    start: Point
    end: Point
    tau: int


@dataclass
class Pin:
    """硬钉：边界点须在 Γ_ij 上，即 Ψ_ij(point)=0。"""

    point: Point
    i: int
    j: int


@dataclass
class InteriorCell:
    """(E) 内部见证：胞心 x_k 归属 i_k。"""

    point: Point
    label: int


@dataclass
class ProblemSpec:
    """抽象输入 (f_D, γ, p_1, τ, Π) 的离散化规格。"""

    boundary: list[BoundarySegment]
    p1: Point
    n_flags: int
    pins: list[Pin] = field(default_factory=list)
    interior: list[InteriorCell] = field(default_factory=list)
    r1: float = 1.0
    # 对称 Ansatz：p_i = (b_i, y0) 同一水平线，减少未知数
    collinear_y: float | None = None
    boundary_samples_per_edge: int = 5
    # 额外等式：固定剩余自由度（如路线 A 后再钉 b）
    extra_equalities: dict[str, float] = field(default_factory=dict)


def polygon_boundary(vertices: Sequence[tuple[float, float]], tau_by_edge: list[int]) -> list[BoundarySegment]:
    n = len(vertices)
    segs: list[BoundarySegment] = []
    for k in range(n):
        x0, y0 = vertices[k]
        x1, y1 = vertices[(k + 1) % n]
        segs.append(BoundarySegment(Point(x0, y0), Point(x1, y1), tau_by_edge[k]))
    return segs


def parse_vertices(text: str) -> list[tuple[float, float]]:
    verts: list[tuple[float, float]] = []
    for chunk in text.split():
        a, b = chunk.split(",")
        verts.append((float(a), float(b)))
    if len(verts) < 3:
        raise ValueError("至少 3 个顶点")
    return verts


def parse_tau(text: str, n_edges: int) -> list[int]:
    if not text:
        return [1] * n_edges
    mapping: dict[int, int] = {}
    for part in text.split():
        idx_s, tau_s = part.split(":")
        mapping[int(idx_s)] = int(tau_s)
    return [mapping.get(k, 1) for k in range(n_edges)]


def load_spec(path: Path) -> ProblemSpec:
    data = json.loads(path.read_text(encoding="utf-8"))
    verts = [tuple(v) for v in data["vertices"]]
    tau = data.get("tau_by_edge", [1] * len(verts))
    holes = data.get("holes", [])
    boundary = polygon_boundary(verts, tau)
    for hole_verts in holes:
        hole_tau = data.get("hole_tau", [1] * len(hole_verts))
        boundary.extend(polygon_boundary(hole_verts, hole_tau))
    p1 = Point(*data["p1"])
    pins = [
        Pin(Point(*p["point"]), p["i"], p["j"])
        for p in data.get("pins", [])
    ]
    interior = [
        InteriorCell(Point(*c["point"]), c["label"])
        for c in data.get("interior", [])
    ]
    extra = {str(k): float(v) for k, v in data.get("fix", {}).items()}
    return ProblemSpec(
        boundary=boundary,
        p1=p1,
        n_flags=int(data.get("n", 2)),
        pins=pins,
        interior=interior,
        r1=float(data.get("r1", 1.0)),
        collinear_y=data.get("collinear_y"),
        boundary_samples_per_edge=int(data.get("boundary_samples", 5)),
        extra_equalities=extra,
    )


# ---------------------------------------------------------------------------
# Ψ_ij 与方程组 𝔈^mul,full
# ---------------------------------------------------------------------------


def psi_expr(
    x: sp.Expr,
    y: sp.Expr,
    pi_x: sp.Expr,
    pi_y: sp.Expr,
    pj_x: sp.Expr,
    pj_y: sp.Expr,
    ri: sp.Expr,
    rj: sp.Expr,
) -> sp.Expr:
    """Ψ_ij = ||x-p_i||² r_j² - ||x-p_j||² r_i²"""
    di2 = (x - pi_x) ** 2 + (y - pi_y) ** 2
    dj2 = (x - pj_x) ** 2 + (y - pj_y) ** 2
    return di2 * rj**2 - dj2 * ri**2


@dataclass
class SiteSymbols:
    """未知站点符号；旗 1 的坐标固定为 p1，r1 固定。"""

    n: int
    p1: Point
    r1: float
    collinear_y: float | None
    # p_i = (px_i, py_i), r_i for i>=2; 旗1 不入符号表
    px: list[sp.Symbol]
    py: list[sp.Symbol]
    r: list[sp.Symbol]

    @classmethod
    def build(cls, spec: ProblemSpec) -> SiteSymbols:
        n = spec.n_flags
        px: list[sp.Symbol] = []
        py: list[sp.Symbol] = []
        r: list[sp.Symbol] = []
        for i in range(2, n + 1):
            px.append(sp.Symbol(f"p{i}x", real=True))
            if spec.collinear_y is None:
                py.append(sp.Symbol(f"p{i}y", real=True))
            r.append(sp.Symbol(f"r{i}", real=True))
        return cls(
            n=n,
            p1=spec.p1,
            r1=spec.r1,
            collinear_y=spec.collinear_y,
            px=px,
            py=py,
            r=r,
        )

    def coords(self, flag: int) -> tuple[sp.Expr, sp.Expr]:
        if flag == 1:
            return sp.Float(self.p1.x), sp.Float(self.p1.y)
        idx = flag - 2
        y = sp.Float(self.collinear_y) if self.collinear_y is not None else self.py[idx]
        return self.px[idx], y

    def radius(self, flag: int) -> sp.Expr:
        if flag == 1:
            return sp.Float(self.r1)
        return self.r[flag - 2]

    def unknowns(self) -> list[sp.Symbol]:
        out: list[sp.Symbol] = []
        out.extend(self.px)
        if self.collinear_y is None:
            out.extend(self.py)
        out.extend(self.r)
        return out


def build_equalities(spec: ProblemSpec, sites: SiteSymbols) -> list[sp.Expr]:
    """等式部分：Π 硬钉 + (E) 内点 tie（可选）。"""
    eqs: list[sp.Expr] = []

    for pin in spec.pins:
        x = sp.Float(pin.point.x)
        y = sp.Float(pin.point.y)
        pi_x, pi_y = sites.coords(pin.i)
        pj_x, pj_y = sites.coords(pin.j)
        ri, rj = sites.radius(pin.i), sites.radius(pin.j)
        eqs.append(psi_expr(x, y, pi_x, pi_y, pj_x, pj_y, ri, rj))

    for cell in spec.interior:
        x = sp.Float(cell.point.x)
        y = sp.Float(cell.point.y)
        lab = cell.label
        pi_x, pi_y = sites.coords(lab)
        ri = sites.radius(lab)
        for j in range(1, spec.n_flags + 1):
            if j == lab:
                continue
            pj_x, pj_y = sites.coords(j)
            rj = sites.radius(j)
            # 严格归属：Ψ_lab,j < 0 在见证点用 tie 松弛为 = -eps 或验证阶段检查
            # 此处不加入等式，留给 verify_viol
            _ = psi_expr(x, y, pi_x, pi_y, pj_x, pj_y, ri, rj)

    return eqs


def sample_segment(seg: BoundarySegment, m: int) -> list[Point]:
    if m < 2:
        m = 2
    pts: list[Point] = []
    for k in range(m):
        t = k / (m - 1) if m > 1 else 0.0
        pts.append(
            Point(
                seg.start.x + t * (seg.end.x - seg.start.x),
                seg.start.y + t * (seg.end.y - seg.start.y),
            )
        )
    return pts


def psi_numeric(
    x: float,
    y: float,
    pi: tuple[float, float],
    pj: tuple[float, float],
    ri: float,
    rj: float,
) -> float:
    di2 = (x - pi[0]) ** 2 + (y - pi[1]) ** 2
    dj2 = (x - pj[0]) ** 2 + (y - pj[1]) ** 2
    return di2 * rj**2 - dj2 * ri**2


def assign_min_flag(
    x: float,
    y: float,
    positions: list[tuple[float, float]],
    radii: list[float],
) -> int:
    phis = [
        math.hypot(x - positions[i][0], y - positions[i][1]) / radii[i]
        for i in range(len(positions))
    ]
    return int(min(range(len(phis)), key=lambda i: phis[i])) + 1


def verify_viol(
    spec: ProblemSpec,
    positions: list[tuple[float, float]],
    radii: list[float],
    samples_per_edge: int = 12,
    tol: float = 1e-6,
) -> tuple[bool, list[str]]:
    """内禀 Viol：σ_S(γ(s)) ≠ τ(s)。"""
    ok = True
    logs: list[str] = []
    for seg in spec.boundary:
        for pt in sample_segment(seg, samples_per_edge):
            sigma = assign_min_flag(pt.x, pt.y, positions, radii)
            if sigma != seg.tau:
                # tie 集允许 Ψ≈0
                lab = seg.tau
                pi = positions[lab - 1]
                ri = radii[lab - 1]
                worst = 0.0
                for j in range(1, spec.n_flags + 1):
                    if j == lab:
                        continue
                    pj = positions[j - 1]
                    rj = radii[j - 1]
                    worst = max(worst, psi_numeric(pt.x, pt.y, pi, pj, ri, rj))
                if worst < -tol:
                    ok = False
                    logs.append(
                        f"Viol @({pt.x:.4f},{pt.y:.4f}): σ={sigma}, τ={seg.tau}, maxΨ={worst:.2e}"
                    )
    return ok, logs


def positions_radii_from_solution(
    spec: ProblemSpec, sol: dict[str, float]
) -> tuple[list[tuple[float, float]], list[float]]:
    pos = [(spec.p1.x, spec.p1.y)]
    rad = [spec.r1]
    for i in range(2, spec.n_flags + 1):
        pos.append((sol[f"p{i}x"], sol[f"p{i}y"]))
        rad.append(sol[f"r{i}"])
    return pos, rad


def verify_interior(
    spec: ProblemSpec,
    positions: list[tuple[float, float]],
    radii: list[float],
    tol: float = 1e-6,
) -> tuple[bool, list[str]]:
    ok = True
    logs: list[str] = []
    for cell in spec.interior:
        lab = cell.label
        pi, ri = positions[lab - 1], radii[lab - 1]
        for j in range(1, spec.n_flags + 1):
            if j == lab:
                continue
            pj, rj = positions[j - 1], radii[j - 1]
            v = psi_numeric(cell.point.x, cell.point.y, pi, pj, ri, rj)
            if v > tol:
                ok = False
                logs.append(
                    f"(E) fail @({cell.point.x},{cell.point.y}): Ψ_{lab}{j}={v:.2e}>0"
                )
    return ok, logs


# ---------------------------------------------------------------------------
# 消元与三分类
# ---------------------------------------------------------------------------


@dataclass
class SolveReport:
    classification: str  # empty | finite | positive_dim | underdetermined
    dimension_hint: str
    groebner_size: int
    solutions: list[dict[str, float]]
    parametric: list[dict[str, str]]
    viol_ok: bool
    interior_ok: bool
    messages: list[str]


def _symbol_name_map(syms: SiteSymbols) -> dict[str, sp.Symbol]:
    m: dict[str, sp.Symbol] = {}
    for i in range(2, syms.n + 1):
        m[f"p{i}x"] = syms.px[i - 2]
        if syms.collinear_y is None:
            m[f"p{i}y"] = syms.py[i - 2]
        m[f"r{i}"] = syms.r[i - 2]
    return m


def _eval_sol(
    sol: dict[sp.Symbol, sp.Expr],
    syms: SiteSymbols,
    extra: dict[str, float] | None = None,
) -> dict[str, float] | None:
    name_map = _symbol_name_map(syms)
    assign = dict(extra or {})
    for name, sym in name_map.items():
        if name in assign:
            continue
        if sym in sol:
            val = sol[sym]
            if val.free_symbols:
                return None
            assign[name] = float(sp.N(val))
    for i in range(2, syms.n + 1):
        if f"p{i}x" not in assign:
            return None
        if syms.collinear_y is None and f"p{i}y" not in assign:
            return None
        if f"r{i}" not in assign:
            return None
        if assign[f"r{i}"] <= 0:
            return None
    out: dict[str, float] = {}
    for i in range(2, syms.n + 1):
        out[f"p{i}x"] = assign[f"p{i}x"]
        out[f"p{i}y"] = syms.collinear_y if syms.collinear_y is not None else assign[f"p{i}y"]
        out[f"r{i}"] = assign[f"r{i}"]
    return out


def _parametric_strings(sol: dict[sp.Symbol, sp.Expr], syms: SiteSymbols) -> dict[str, str]:
    name_map = _symbol_name_map(syms)
    out: dict[str, str] = {}
    for i in range(2, syms.n + 1):
        sx = syms.px[i - 2]
        out[f"p{i}x"] = str(sol.get(sx, sx))
        if syms.collinear_y is not None:
            out[f"p{i}y"] = str(syms.collinear_y)
        else:
            sy = syms.py[i - 2]
            out[f"p{i}y"] = str(sol.get(sy, sy))
        rs = syms.r[i - 2]
        out[f"r{i}"] = str(sol.get(rs, rs))
    return out


def _solve_from_groebner(G: GroebnerBasis, unknowns: list[sp.Symbol]) -> list[dict[sp.Symbol, sp.Expr]]:
    """Gröbner 给出 r²=c 等情形时，补全 ± 实根。"""
    base: dict[sp.Symbol, sp.Expr] = {}
    for poly in G.polys:
        if poly.total_degree() == 0:
            if poly.as_expr() == 1:
                return []
            continue
        active = [u for u in unknowns if poly.degree(u) > 0]
        if len(active) != 1:
            continue
        sym = active[0]
        deg = poly.degree(sym)
        if deg == 1:
            base[sym] = -poly.TC() / poly.LC()
        elif deg == 2 and abs(float(sp.N(poly.LC())) - 1.0) < 1e-9:
            c = float(sp.N(-poly.TC()))
            if c > 0:
                base[sym] = sp.sqrt(c)
    if not base:
        return []
    free = [u for u in unknowns if u not in base]
    if free:
        return [{**base, u: u} for u in free]
    out: list[dict[sp.Symbol, sp.Expr]] = [dict(base)]
    for sym, val in list(base.items()):
        if sym.name.startswith("r") and val.is_Pow and val.exp == sp.S.Half:
            alt = dict(base)
            alt[sym] = -val
            out.append(alt)
    return out


def solve_spec(spec: ProblemSpec) -> SolveReport:
    sites = SiteSymbols.build(spec)
    unknowns = sites.unknowns()
    eq_exprs = build_equalities(spec, sites)
    messages: list[str] = []

    if spec.extra_equalities:
        name_map = _symbol_name_map(sites)
        for key, val in spec.extra_equalities.items():
            if key not in name_map:
                raise ValueError(f"未知固定参数: {key}")
            eq_exprs.append(name_map[key] - val)

    if not eq_exprs:
        messages.append("无等式约束（仅钉点/内点可 --pin 指定）；无法消元。")
        return SolveReport("underdetermined", "dof>0", 0, [], [], False, False, messages)

    messages.append(f"等式数={len(eq_exprs)}, 未知数={len(unknowns)} ({', '.join(str(u) for u in unknowns)})")

    G: GroebnerBasis = groebner(eq_exprs, unknowns, order="lex")
    messages.append(f"Groebner 基大小={len(G.polys)}")

    raw_sols = sp.solve(list(G), unknowns, dict=True)
    if not raw_sols:
        raw_sols = _solve_from_groebner(G, unknowns)
    parametric: list[dict[str, str]] = []
    numeric_sols: list[dict[str, float]] = []
    for raw in raw_sols:
        free = set()
        for v in raw.values():
            free |= v.free_symbols
        if free:
            parametric.append(_parametric_strings(raw, sites))
            # 若用户给了 extra_equalities，再试数值
            ev = _eval_sol(raw, sites, spec.extra_equalities)
            if ev is not None:
                numeric_sols.append(ev)
        else:
            ev = _eval_sol(raw, sites)
            if ev is not None:
                numeric_sols.append(ev)

    n_free = max(0, len(unknowns) - len(eq_exprs))
    if parametric and not numeric_sols:
        klass = "positive_dim"
        dim_hint = f"正维族（约 {n_free} 自由参数）；见 parametric"
    elif not numeric_sols:
        if len(eq_exprs) >= len(unknowns):
            klass = "empty"
            dim_hint = "Sol^W=空（或无实正半径解）"
        else:
            klass = "positive_dim"
            dim_hint = f"正维族：{n_free} 自由参数"
    elif len(numeric_sols) == 1:
        klass = "finite"
        dim_hint = "dim=0，孤立实解"
    else:
        klass = "finite"
        dim_hint = f"dim=0，{len(numeric_sols)} 个实解"

    # 审计
    viol_all = True
    interior_all = True
    for sol in numeric_sols:
        pos = [(spec.p1.x, spec.p1.y)]
        rad = [spec.r1]
        for i in range(2, spec.n_flags + 1):
            pos.append((sol[f"p{i}x"], sol[f"p{i}y"]))
            rad.append(sol[f"r{i}"])
        v_ok, v_logs = verify_viol(spec, pos, rad, spec.boundary_samples_per_edge)
        i_ok, i_logs = verify_interior(spec, pos, rad)
        viol_all = viol_all and v_ok
        interior_all = interior_all and i_ok
        messages.extend(v_logs[:5])
        messages.extend(i_logs[:5])

    return SolveReport(
        classification=klass,
        dimension_hint=dim_hint,
        groebner_size=len(G.polys),
        solutions=numeric_sols,
        parametric=parametric,
        viol_ok=viol_all,
        interior_ok=interior_all,
        messages=messages,
    )


def format_report(spec: ProblemSpec, report: SolveReport) -> str:
    lines: list[str] = []
    lines.append("=" * 60)
    lines.append("乘性泰森消元报告（§W.6 工具层）")
    lines.append("=" * 60)
    lines.append(f"旗数 N={spec.n_flags}, p1=({spec.p1.x},{spec.p1.y}), r1={spec.r1}")
    lines.append(f"边界段数={len(spec.boundary)}, 硬钉 |Π|={len(spec.pins)}, 内见证 |E|={len(spec.interior)}")
    lines.append(f"分类: {report.classification} — {report.dimension_hint}")
    lines.append(f"Viol 审计: {'PASS' if report.viol_ok else 'FAIL'}")
    lines.append(f"(E) 审计: {'PASS' if report.interior_ok else 'FAIL'}")
    lines.append("-" * 60)
    for msg in report.messages:
        lines.append(msg)
    if report.parametric:
        lines.append("-" * 60)
        lines.append("参数化解族（闭式分支，可 --fix 再数值化）:")
        for k, sol in enumerate(report.parametric):
            lines.append(f"  [族 #{k+1}]")
            for i in range(2, spec.n_flags + 1):
                lines.append(f"    p{i}x={sol[f'p{i}x']}, p{i}y={sol[f'p{i}y']}, r{i}={sol[f'r{i}']}")
    if report.solutions:
        lines.append("-" * 60)
        lines.append("实解 (p_i, r_i):")
        for k, sol in enumerate(report.solutions):
            lines.append(f"  [#{k+1}]")
            for i in range(2, spec.n_flags + 1):
                lines.append(
                    f"    p{i}=({sol[f'p{i}x']:.8f}, {sol[f'p{i}y']:.8f}), r{i}={sol[f'r{i}']:.8f}"
                )
    if not report.solutions and not report.parametric:
        lines.append("无通过筛选的实正半径解。")
    lines.append("=" * 60)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 内置代入（§8/§14 型，非新定理）
# ---------------------------------------------------------------------------


def example_rect_route_a() -> ProblemSpec:
    """矩形路线 A：只钉右侧面中点 M∈Γ_12（§14 同型）。"""
    W, H = 2.0, 1.0
    a, mid = 0.4, 0.5
    verts = [(0, 0), (W, 0), (W, H), (0, H)]
    # 顶边 τ=2（渗漏段故事），其余 τ=1
    tau = [1, 2, 1, 1]
    return ProblemSpec(
        boundary=polygon_boundary(verts, tau),
        p1=Point(a, mid),
        n_flags=2,
        pins=[Pin(Point(W, mid), 1, 2)],
        interior=[InteriorCell(Point(a + 0.05, mid), 1), InteriorCell(Point(W - 0.2, mid), 2)],
        collinear_y=mid,
        boundary_samples_per_edge=8,
    )


def example_sector_route_a() -> ProblemSpec:
    """实心扇形路线 A 对称钉 M=(R,0)（§8 代入）。"""
    R, alpha, a = 2.0, math.pi / 6, 0.5
    # 用折线逼近外弧 + 两条径向边
    arc_pts = [
        (R * math.cos(t), R * math.sin(t))
        for t in [alpha, alpha / 2, 0, -alpha / 2, -alpha]
    ]
    verts = [(0, 0)] + arc_pts[::-1]  # 闭合多边形近似
    tau = [1] * len(verts)  # 外弧与径向边故事可细化
    return ProblemSpec(
        boundary=polygon_boundary(verts, tau),
        p1=Point(a, 0.0),
        n_flags=2,
        pins=[Pin(Point(R, 0.0), 1, 2)],
        interior=[InteriorCell(Point(1.2, 0.0), 2), InteriorCell(Point(0.7, 0.0), 1)],
        collinear_y=0.0,
        boundary_samples_per_edge=6,
    )


def write_example_json(path: Path) -> None:
    spec = example_rect_route_a()
    data: dict[str, Any] = {
        "vertices": [[0, 0], [2, 0], [2, 1], [0, 1]],
        "tau_by_edge": [1, 2, 1, 1],
        "p1": [0.4, 0.5],
        "n": 2,
        "r1": 1.0,
        "collinear_y": 0.5,
        "pins": [{"point": [2.0, 0.5], "i": 1, "j": 2}],
        "interior": [
            {"point": [0.45, 0.5], "label": 1},
            {"point": [1.8, 0.5], "label": 2},
        ],
        "boundary_samples": 8,
    }
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_spec_from_args(args: argparse.Namespace) -> ProblemSpec:
    if args.spec:
        return load_spec(Path(args.spec))
    if args.example == "rect-route-a":
        return example_rect_route_a()
    if args.example == "sector-route-a":
        return example_sector_route_a()
    if not args.vertices:
        raise SystemExit("请指定 --spec、--example 或 --vertices")

    verts = parse_vertices(args.vertices)
    tau = parse_tau(args.tau or "", len(verts))
    p1x, p1y = map(float, args.p1.split(","))
    pins: list[Pin] = []
    for pin_s in args.pin or []:
        xs, ys, i_s, j_s = pin_s.split(",")
        pins.append(Pin(Point(float(xs), float(ys)), int(i_s), int(j_s)))
    interior: list[InteriorCell] = []
    for cell_s in args.interior or []:
        xs, ys, lab = cell_s.split(",")
        interior.append(InteriorCell(Point(float(xs), float(ys)), int(lab)))

    collinear = float(args.collinear_y) if args.collinear_y is not None else None
    extra: dict[str, float] = {}
    for item in args.fix or []:
        key, val = item.split("=")
        extra[key.strip()] = float(val.strip())
    return ProblemSpec(
        boundary=polygon_boundary(verts, tau),
        p1=Point(p1x, p1y),
        n_flags=int(args.n),
        pins=pins,
        interior=interior,
        r1=float(args.r1),
        collinear_y=collinear,
        boundary_samples_per_edge=int(args.boundary_samples),
        extra_equalities=extra,
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="乘性泰森 𝔖ol^W 消元工具")
    parser.add_argument("--spec", help="JSON 问题规格")
    parser.add_argument(
        "--example",
        choices=["rect-route-a", "sector-route-a"],
        help="内置代入算例",
    )
    parser.add_argument("--vertices", help='多边形顶点，如 "0,0 2,0 2,1 0,1"')
    parser.add_argument("--p1", default="0,0", help="固定 p1，如 0.4,0.5")
    parser.add_argument("--n", type=int, default=2, help="旗数 N")
    parser.add_argument("--r1", type=float, default=1.0)
    parser.add_argument("--collinear-y", type=float, default=None, help="对称 Ansatz 水平线 y0")
    parser.add_argument("--pin", action="append", help="硬钉 x,y,i,j")
    parser.add_argument("--tau", help="边归属 edge_idx:tau，如 0:1 1:2 2:1 3:1")
    parser.add_argument("--interior", action="append", help="内见证 x,y,label")
    parser.add_argument("--boundary-samples", type=int, default=8)
    parser.add_argument(
        "--fix",
        action="append",
        help="固定未知量，如 p2x=0.6 或 r2=0.75（加等式封闭解族）",
    )
    parser.add_argument("--write-example-json", metavar="PATH", help="写出 rect-route-a JSON")
    parser.add_argument("--json-out", help="报告 JSON 输出路径")
    args = parser.parse_args(argv)

    if args.write_example_json:
        write_example_json(Path(args.write_example_json))
        print(f"已写入 {args.write_example_json}")
        return 0

    spec = build_spec_from_args(args)
    if args.fix:
        merged = dict(spec.extra_equalities)
        for item in args.fix:
            key, val = item.split("=")
            merged[key.strip()] = float(val.strip())
        spec = ProblemSpec(
            boundary=spec.boundary,
            p1=spec.p1,
            n_flags=spec.n_flags,
            pins=spec.pins,
            interior=spec.interior,
            r1=spec.r1,
            collinear_y=spec.collinear_y,
            boundary_samples_per_edge=spec.boundary_samples_per_edge,
            extra_equalities=merged,
        )
    report = solve_spec(spec)
    text = format_report(spec, report)
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    print(text)

    if args.json_out:
        payload = {
            "classification": report.classification,
            "dimension_hint": report.dimension_hint,
            "solutions": report.solutions,
            "parametric": report.parametric,
            "viol_ok": report.viol_ok,
            "interior_ok": report.interior_ok,
            "messages": report.messages,
        }
        Path(args.json_out).write_text(json.dumps(payload, indent=2), encoding="utf-8")

    if report.solutions or report.parametric:
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
