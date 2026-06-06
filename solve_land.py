#!/usr/bin/env python3
"""地皮划分 —— 原问题实用封口版。

问题定义（用户规格）：
  输入：
    ① 地皮边界（+ 可选孔洞）
    ② 第一面旗 p1
    ③ 边界归属 τ（每条边归哪面旗）— 显式输入
  过程：
    ④ 自动判定 n_min（定理 M★：m_eff=0→2，m_eff≥1→3 秩3）
  输出：
    · 无解
    · 有限多解 → 列出全部可行解
    · 无穷多解 → 给出解空间（参数化闭式 + 自由参说明）

用法:
  python solve_land.py --vertices 0,0 2,0 2,1 0,1 --p1 0.4,0.5 --tau 0:1 1:2 2:1 3:1
  python solve_land.py --vertices 0,0 3,0 3,1.5 1,1.5 1,3 0,3 --p1 1.25,1.125 --tau 0:1 1:2 2:1 3:1 4:1 5:1 -o out.png
  python solve_land.py --spec examples/rect_route_a.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import sympy as sp

from mul_tyson_solve import (
    BoundarySegment,
    InteriorCell,
    Pin,
    Point,
    ProblemSpec,
    SolveReport,
    load_spec,
    parse_tau,
    parse_vertices,
    polygon_boundary,
    positions_radii_from_solution,
    sample_segment,
    solve_spec,
    verify_interior,
    verify_viol,
)
from mul_tyson_viz import (
    DomainInput,
    combined_score,
    compute_region_grid,
    domain_from_spec_json,
    edge_midpoint,
    find_interior_point,
    instantiate_parametric,
    plot_tyson,
    resolve_solution,
)


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------


@dataclass
class SolutionSpace:
    kind: str  # empty | unique | finite | positive_dim
    n_min: int
    m_eff: int
    rank3: bool
    finite_solutions: list[dict[str, float]] = field(default_factory=list)
    parametric_families: list[dict[str, str]] = field(default_factory=list)
    free_params: list[list[str]] = field(default_factory=list)
    viol_pass_count: int = 0


@dataclass
class LandResult:
    spec: ProblemSpec
    domain: DomainInput
    space: SolutionSpace
    lines: list[str] = field(default_factory=list)

    def text(self) -> str:
        return "\n".join(self.lines)


# ---------------------------------------------------------------------------
# n_min 与 m_eff（定理 M★ 工具层）
# ---------------------------------------------------------------------------


def _chi12(x: float, y: float, p1: Point, p2: tuple[float, float], lam: float) -> float:
    return math.hypot(x - p2[0], y - p2[1]) / lam - math.hypot(x - p1.x, y - p1.y)


def audit_m_eff(
    outer: Sequence[tuple[float, float]],
    p1: Point,
    p2: tuple[float, float],
    lam: float,
    samples: int = 20,
) -> tuple[int, list[int]]:
    leaky: list[int] = []
    n = len(outer)
    for k in range(n):
        seg = BoundarySegment(
            Point(*outer[k]), Point(*outer[(k + 1) % n]), 1
        )
        chi_min = min(
            _chi12(pt.x, pt.y, p1, p2, lam) for pt in sample_segment(seg, samples)
        )
        if chi_min < -1e-5:
            leaky.append(k)
    return len(leaky), leaky


def n_min_from_m_eff(m_eff: int) -> int:
    return 2 if m_eff == 0 else 3


# ---------------------------------------------------------------------------
# 建题：τ 输入 + 自动钉点
# ---------------------------------------------------------------------------


def pin_on_tau_edge(
    vertices: Sequence[tuple[float, float]],
    tau: list[int],
    flag_i: int = 1,
    flag_j: int = 2,
) -> Pin:
    """在第一条 τ=flag_j 的边中点钉 Γ_ij。"""
    for k, t in enumerate(tau):
        if t == flag_j:
            a = vertices[k]
            b = vertices[(k + 1) % len(vertices)]
            return Pin(edge_midpoint(a, b), flag_i, flag_j)
    # 回退：最长边中点
    k = 0
    a = vertices[k]
    b = vertices[(k + 1) % len(vertices)]
    return Pin(edge_midpoint(a, b), flag_i, flag_j)


def interior_witnesses(p1: Point, pin: Pin) -> list[InteriorCell]:
    mx, my = pin.point.x, pin.point.y
    return [
        InteriorCell(Point(p1.x + 0.12 * (p1.x - mx), p1.y + 0.12 * (p1.y - my)), 1),
        InteriorCell(Point(p1.x + 0.5 * (mx - p1.x), p1.y + 0.5 * (my - p1.y)), 2),
    ]


def build_spec(
    domain: DomainInput,
    p1: Point,
    tau: list[int],
    n_flags: int = 2,
    collinear_y: float | None = None,
) -> ProblemSpec:
    boundary = polygon_boundary(domain.vertices, tau)
    for hole in domain.holes:
        boundary.extend(polygon_boundary(hole, [1] * len(hole)))
    pin = pin_on_tau_edge(domain.vertices, tau)
    if collinear_y is None:
        xs = [v[0] for v in domain.vertices]
        ys = [v[1] for v in domain.vertices]
        if max(xs) - min(xs) > 1.2 * max(max(ys) - min(ys), 1e-6):
            collinear_y = p1.y
    return ProblemSpec(
        boundary=boundary,
        p1=p1,
        n_flags=n_flags,
        pins=[pin],
        interior=interior_witnesses(p1, pin),
        collinear_y=collinear_y,
        boundary_samples_per_edge=12,
    )


# ---------------------------------------------------------------------------
# 解空间：有限列举 / 无穷参数化
# ---------------------------------------------------------------------------


def _positions_radii(spec: ProblemSpec, sol: dict[str, float], rank3: bool):
    if rank3:
        p2 = (sol["p2x"], sol["p2y"])
        r2 = sol["r2"]
        return [(spec.p1.x, spec.p1.y), p2, p2], [spec.r1, r2, r2]
    return positions_radii_from_solution(spec, sol)


def verify_solution(
    spec: ProblemSpec, sol: dict[str, float], rank3: bool
) -> tuple[bool, bool]:
    pos, rad = _positions_radii(spec, sol, rank3)
    v_ok, _ = verify_viol(spec, pos, rad, spec.boundary_samples_per_edge)
    i_ok, _ = verify_interior(spec, pos, rad)
    return v_ok, i_ok


def _dedupe_solutions(sols: list[dict[str, float]], tol: float = 1e-4) -> list[dict[str, float]]:
    out: list[dict[str, float]] = []
    for s in sols:
        dup = False
        for t in out:
            if all(abs(s.get(k, 0) - t.get(k, 0)) < tol for k in s):
                dup = True
                break
        if not dup:
            out.append(s)
    return out


def enumerate_groebner_finite(report: SolveReport) -> list[dict[str, float]]:
    return _dedupe_solutions(list(report.solutions))


def free_symbols_in_family(fam: dict[str, str]) -> list[str]:
    free: set[str] = set()
    for text in fam.values():
        free |= {str(s) for s in sp.sympify(text).free_symbols}
    return sorted(free)


def sample_parametric_feasible(
    spec: ProblemSpec,
    fam: dict[str, str],
    rank3: bool,
    n_grid: int = 100,
) -> list[dict[str, float]]:
    """从 1 维族密采样，收集 Viol PASS 的孤立点（展示多解切片）。"""
    free = free_symbols_in_family(fam)
    if not free:
        inst = instantiate_parametric(fam, {})
        return [inst] if inst and verify_solution(spec, inst, rank3)[0] else []
    if len(free) > 2:
        return []
    found: list[dict[str, float]] = []
    grids: list[np.ndarray]
    if len(free) == 1:
        name = free[0]
        lo, hi = (0.08, 3.5) if name.startswith("r") else (-2.0, 5.0)
        grids = [np.linspace(lo, hi, n_grid)]
        keys = [name]
    else:
        keys = free[:2]
        grids = [np.linspace(0.08, 3.5, 40), np.linspace(-1.0, 4.0, 40)]

    if len(free) == 1:
        for val in grids[0]:
            inst = instantiate_parametric(fam, {keys[0]: float(val)})
            if inst and verify_solution(spec, inst, rank3)[0]:
                found.append(inst)
    else:
        for v1 in grids[0]:
            for v2 in grids[1]:
                inst = instantiate_parametric(fam, {keys[0]: float(v1), keys[1]: float(v2)})
                if inst and verify_solution(spec, inst, rank3)[0]:
                    found.append(inst)
    return _dedupe_solutions(found)


def classify_solution_space(
    spec: ProblemSpec,
    report: SolveReport,
    domain: DomainInput,
    rank3: bool,
) -> SolutionSpace:
    finite = enumerate_groebner_finite(report)
    feasible_finite: list[dict[str, float]] = []
    for sol in finite:
        if verify_solution(spec, sol, rank3)[0]:
            feasible_finite.append(sol)

    families = list(report.parametric)
    free_lists = [free_symbols_in_family(f) for f in families]

    # 从参数族采样补充（Gröbner 有时只给族不给孤立点）
    for fam in families:
        for sol in sample_parametric_feasible(spec, fam, rank3, n_grid=80):
            if sol not in feasible_finite:
                feasible_finite.append(sol)
    feasible_finite = _dedupe_solutions(feasible_finite)

    # m_eff 用第一个可行或首个有限解估计
    m_eff, _ = 0, []
    if feasible_finite:
        s = feasible_finite[0]
        m_eff, _ = audit_m_eff(
            domain.vertices,
            spec.p1,
            (s["p2x"], s["p2y"]),
            s["r2"],
        )
    elif finite:
        s = finite[0]
        m_eff, _ = audit_m_eff(
            domain.vertices, spec.p1, (s["p2x"], s["p2y"]), s["r2"]
        )

    n_min = n_min_from_m_eff(m_eff)

    if not feasible_finite and not families:
        kind = "empty"
    elif families and not feasible_finite:
        kind = "positive_dim"
    elif families and feasible_finite:
        kind = "positive_dim"  # 族 + 部分孤立点仍属无穷
    elif len(feasible_finite) == 1:
        kind = "unique"
    else:
        kind = "finite"

    # 若只有参数族、无孤立可行点
    if kind == "positive_dim" and not families:
        kind = "empty"

    viol_n = len(feasible_finite)
    if kind == "positive_dim":
        for fam in families:
            viol_n += len(sample_parametric_feasible(spec, fam, rank3, n_grid=30))

    return SolutionSpace(
        kind=kind,
        n_min=n_min,
        m_eff=m_eff,
        rank3=rank3,
        finite_solutions=feasible_finite,
        parametric_families=families,
        free_params=free_lists,
        viol_pass_count=viol_n,
    )


def format_solution(sol: dict[str, float], rank3: bool, idx: int) -> str:
    p2x, p2y, r2 = sol["p2x"], sol["p2y"], sol["r2"]
    lines = [f"  解 #{idx}: p2=({p2x:.6f}, {p2y:.6f}), r2={r2:.6f}"]
    if rank3:
        lines.append(f"         p3=p2, r3=r2（秩3）")
    return "\n".join(lines)


def format_parametric_family(fam: dict[str, str], free: list[str], idx: int) -> str:
    lines = [f"  族 #{idx}（自由参: {', '.join(free) if free else '无'}）:"]
    for i in range(2, 6):
        if f"p{i}x" in fam:
            lines.append(f"    p{i}x = {fam[f'p{i}x']}")
            lines.append(f"    p{i}y = {fam.get(f'p{i}y', '?')}")
            lines.append(f"    r{i}  = {fam.get(f'r{i}', '?')}")
    if not any(f"p{i}x" in fam for i in range(2, 6)):
        for k, v in fam.items():
            lines.append(f"    {k} = {v}")
    return "\n".join(lines)


def solve_land_problem(
    domain: DomainInput,
    p1: Point,
    tau: list[int],
) -> LandResult:
    spec = build_spec(domain, p1, tau)
    report = solve_spec(spec)

    # 先按 N=2 消元；m_eff 判定 n_min；m_eff≥1 则秩3 验证
    probe_sol = report.solutions[0] if report.solutions else None
    if probe_sol is None and report.parametric:
        inst = instantiate_parametric(report.parametric[0], {})
        if inst:
            probe_sol = inst
    m_eff = 0
    if probe_sol:
        m_eff, _ = audit_m_eff(
            domain.vertices,
            p1,
            (probe_sol["p2x"], probe_sol["p2y"]),
            probe_sol["r2"],
        )
    rank3 = m_eff >= 1
    n_min = n_min_from_m_eff(m_eff)

    space = classify_solution_space(spec, report, domain, rank3)

    lines: list[str] = []
    lines.append("=" * 58)
    lines.append("【地皮划分 — 原问题封口版】")
    lines.append("=" * 58)
    lines.append("输入:")
    lines.append(f"  地皮：{len(domain.vertices)} 顶点" + (f"，{len(domain.holes)} 孔" if domain.holes else ""))
    lines.append(f"  p1 = ({p1.x:.4f}, {p1.y:.4f})")
    lines.append(f"  τ  = {tau}")
    lines.append("")
    lines.append("n_min 判定（定理 M★，秩3 路线）:")
    lines.append(f"  m_eff = {space.m_eff}  →  n_min = {space.n_min}" + ("（需秩3 第三旗）" if rank3 else ""))
    lines.append("")

    if space.kind == "empty":
        lines.append("结论：无解（𝔖ol^W = ∅ 或无 Viol PASS 实解）")
        if report.parametric:
            lines.append("")
            lines.append("代数解族存在，但边界故事 τ 下无可行实解：")
            for i, fam in enumerate(space.parametric_families):
                lines.append(format_parametric_family(fam, space.free_params[i] if i < len(space.free_params) else [], i + 1))
    elif space.kind == "unique":
        lines.append("结论：唯一解（有限且 Viol PASS 仅 1 个）")
        lines.append("")
        lines.append("解列表:")
        for i, sol in enumerate(space.finite_solutions, 1):
            lines.append(format_solution(sol, rank3, i))
    elif space.kind == "finite":
        lines.append(f"结论：多解（有限 {len(space.finite_solutions)} 个，全部列出）")
        lines.append("")
        lines.append("解列表:")
        for i, sol in enumerate(space.finite_solutions, 1):
            lines.append(format_solution(sol, rank3, i))
    else:  # positive_dim
        lines.append("结论：多解（无穷 — 正维解空间）")
        lines.append("")
        lines.append("解空间（参数化闭式）:")
        for i, fam in enumerate(space.parametric_families):
            free = space.free_params[i] if i < len(space.free_params) else []
            lines.append(format_parametric_family(fam, free, i + 1))
        if space.finite_solutions:
            lines.append("")
            lines.append(
                f"注：解族上 Viol PASS（连续 1 维族，不逐点列举；"
                f"采样检出 {len(space.finite_solutions)} 个代表点）。"
            )
            lines.append("  代表点示例：")
            for i, sol in enumerate(space.finite_solutions[:3], 1):
                lines.append(format_solution(sol, rank3, i))

    lines.append("=" * 58)
    return LandResult(spec=spec, domain=domain, space=space, lines=lines)


def pick_plot_solution(
    spec: ProblemSpec, report: SolveReport, space: SolutionSpace
) -> dict[str, float] | None:
    sol, _ = resolve_solution(spec, report)
    if sol is not None:
        return sol
    if space.finite_solutions:
        return min(space.finite_solutions, key=lambda s: combined_score(spec, s))
    return None


def plot_result(
    domain: DomainInput,
    spec: ProblemSpec,
    sol: dict[str, float],
    rank3: bool,
    out: Path,
) -> None:
    if rank3:
        ext = ProblemSpec(
            boundary=spec.boundary,
            p1=spec.p1,
            n_flags=3,
            pins=spec.pins,
            interior=spec.interior,
            r1=spec.r1,
            collinear_y=spec.collinear_y,
            boundary_samples_per_edge=spec.boundary_samples_per_edge,
        )
        ext_sol = {
            "p2x": sol["p2x"],
            "p2y": sol["p2y"],
            "r2": sol["r2"],
            "p3x": sol["p2x"],
            "p3y": sol["p2y"],
            "r3": sol["r2"],
        }
        plot_tyson(domain, ext, ext_sol, out)
    else:
        plot_tyson(domain, spec, sol, out)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="地皮划分原问题求解器")
    parser.add_argument("--vertices", help='外边界 "x,y x,y ..."')
    parser.add_argument("--holes", action="append")
    parser.add_argument("--spec", help="JSON（含 vertices, tau, p1）")
    parser.add_argument("--p1", help="第一面旗")
    parser.add_argument("--tau", required=False, help="边归属 0:1 1:2 ...（--spec 可省略）")
    parser.add_argument("-o", "--out", help="输出 PNG（取第一个可行解）")
    parser.add_argument("--json-out", help="输出 JSON 报告")
    args = parser.parse_args(argv)

    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    if args.spec:
        path = Path(args.spec)
        domain = domain_from_spec_json(path)
        raw = json.loads(path.read_text(encoding="utf-8"))
        p1 = Point(*raw["p1"])
        tau = raw.get("tau_by_edge", [1] * len(domain.vertices))
        if args.p1:
            p1 = Point(*map(float, args.p1.split(",")))
        if args.tau:
            tau = parse_tau(args.tau, len(domain.vertices))
    elif args.vertices:
        verts = parse_vertices(args.vertices)
        holes = [parse_vertices(h) for h in (args.holes or [])]
        domain = DomainInput(vertices=verts, holes=holes)
        if not args.tau:
            parser.error("需要 --tau（边界归属）")
        tau = parse_tau(args.tau, len(verts))
        p1 = find_interior_point(verts, holes)
        if args.p1:
            p1 = Point(*map(float, args.p1.split(",")))
    else:
        parser.error("需要 --vertices 或 --spec")

    result = solve_land_problem(domain, p1, tau)
    print(result.text())

    if args.json_out:
        payload: dict[str, Any] = {
            "kind": result.space.kind,
            "n_min": result.space.n_min,
            "m_eff": result.space.m_eff,
            "rank3": result.space.rank3,
            "tau": tau,
            "p1": [p1.x, p1.y],
            "finite_solutions": result.space.finite_solutions,
            "parametric_families": result.space.parametric_families,
            "free_params": result.space.free_params,
        }
        Path(args.json_out).write_text(
            json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
        )

    if args.out:
        report = solve_spec(result.spec)
        sol = pick_plot_solution(result.spec, report, result.space)
        if sol:
            plot_result(domain, result.spec, sol, result.space.rank3, Path(args.out))
            print(f"\n已保存图（代表解）: {Path(args.out).resolve()}")
        else:
            print("\n无可行解可绘图。")
            return 1

    return 0 if result.space.kind != "empty" else 1


if __name__ == "__main__":
    raise SystemExit(main())
