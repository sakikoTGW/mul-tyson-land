#!/usr/bin/env python3
"""任意域乘性泰森：自动建题 → 消元 → 可视化。

用法:
  python mul_tyson_viz.py --vertices 0,0 3,0 3,1.5 1,1.5 1,3 0,3
  python mul_tyson_viz.py --vertices 0,0 2,0 2,1 0,1 --p1 0.4,0.5 --out out.png
  python mul_tyson_viz.py --spec examples/rect_route_a.json --out rect.png
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection

from mul_tyson_solve import (
    BoundarySegment,
    InteriorCell,
    Pin,
    Point,
    ProblemSpec,
    SolveReport,
    assign_min_flag,
    format_report,
    load_spec,
    parse_tau,
    parse_vertices,
    polygon_boundary,
    positions_radii_from_solution,
    solve_spec,
    verify_interior,
    verify_viol,
)


# ---------------------------------------------------------------------------
# 几何：点在多边形内（含孔洞）
# ---------------------------------------------------------------------------


def point_in_ring(x: float, y: float, ring: Sequence[tuple[float, float]]) -> bool:
    inside = False
    n = len(ring)
    for i in range(n):
        x0, y0 = ring[i]
        x1, y1 = ring[(i + 1) % n]
        if ((y0 > y) != (y1 > y)) and (
            x < (x1 - x0) * (y - y0) / (y1 - y0 + 1e-15) + x0
        ):
            inside = not inside
    return inside


def point_in_domain(
    x: float,
    y: float,
    outer: Sequence[tuple[float, float]],
    holes: Sequence[Sequence[tuple[float, float]]] = (),
) -> bool:
    if not point_in_ring(x, y, outer):
        return False
    return all(not point_in_ring(x, y, h) for h in holes)


def polygon_centroid(ring: Sequence[tuple[float, float]]) -> Point:
    a = 0.0
    cx = 0.0
    cy = 0.0
    n = len(ring)
    for i in range(n):
        x0, y0 = ring[i]
        x1, y1 = ring[(i + 1) % n]
        cross = x0 * y1 - x1 * y0
        a += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    if abs(a) < 1e-12:
        xs = [p[0] for p in ring]
        ys = [p[1] for p in ring]
        return Point(sum(xs) / n, sum(ys) / n)
    a *= 0.5
    return Point(cx / (6 * a), cy / (6 * a))


def bbox(ring: Sequence[tuple[float, float]], pad: float = 0.05) -> tuple[float, float, float, float]:
    xs = [p[0] for p in ring]
    ys = [p[1] for p in ring]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    w, h = xmax - xmin, ymax - ymin
    return xmin - pad * w, xmax + pad * w, ymin - pad * h, ymax + pad * h


def find_interior_point(
    outer: Sequence[tuple[float, float]],
    holes: Sequence[Sequence[tuple[float, float]]] = (),
) -> Point:
    c = polygon_centroid(outer)
    if point_in_domain(c.x, c.y, outer, holes):
        return c
    x0, y0 = outer[0]
    for t in np.linspace(0.05, 0.95, 40):
        px = x0 + t * (c.x - x0)
        py = y0 + t * (c.y - y0)
        if point_in_domain(px, py, outer, holes):
            return Point(float(px), float(py))
    # 栅格搜索
    xmin, xmax, ymin, ymax = bbox(outer, 0.0)
    for _ in range(5):
        nx, ny = 30, 30
        for ix in range(1, nx - 1):
            for iy in range(1, ny - 1):
                px = xmin + (xmax - xmin) * ix / (nx - 1)
                py = ymin + (ymax - ymin) * iy / (ny - 1)
                if point_in_domain(px, py, outer, holes):
                    return Point(float(px), float(py))
    raise ValueError("无法在域内找到内点")


def edge_midpoint(a: tuple[float, float], b: tuple[float, float]) -> Point:
    return Point(0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]))


def edge_distance_to_point(
    a: tuple[float, float], b: tuple[float, float], p: Point
) -> float:
    m = edge_midpoint(a, b)
    return math.hypot(m.x - p.x, m.y - p.y)


# ---------------------------------------------------------------------------
# 自动建题（路线 A 默认故事）
# ---------------------------------------------------------------------------


@dataclass
class DomainInput:
    vertices: list[tuple[float, float]]
    holes: list[list[tuple[float, float]]] = field(default_factory=list)


def auto_story_tau(
    vertices: Sequence[tuple[float, float]], p1: Point, n_flags: int = 2
) -> list[int]:
    """最远外特征边归旗 2，其余归旗 1（N=2 默认故事）。"""
    if n_flags != 2:
        return [1] * len(vertices)
    dists = []
    for i in range(len(vertices)):
        a = vertices[i]
        b = vertices[(i + 1) % len(vertices)]
        dists.append(edge_distance_to_point(a, b, p1))
    far_idx = int(max(range(len(dists)), key=lambda i: dists[i]))
    tau = [1] * len(vertices)
    tau[far_idx] = 2
    return tau


def auto_pin_route_a(
    vertices: Sequence[tuple[float, float]], tau: list[int], i: int = 1, j: int = 2
) -> Pin:
    """在 τ=2 的边中点钉 Γ_ij（若无 τ=2 则取最远边）。"""
    candidates = [k for k, t in enumerate(tau) if t == j]
    if not candidates:
        candidates = [0]
    k = candidates[0]
    a = vertices[k]
    b = vertices[(k + 1) % len(vertices)]
    return Pin(edge_midpoint(a, b), i, j)


def auto_interior_witnesses(
    p1: Point, pin: Pin, n_flags: int = 2
) -> list[InteriorCell]:
    if n_flags < 2:
        return []
    mx, my = pin.point.x, pin.point.y
    w1 = Point(p1.x + 0.15 * (p1.x - mx), p1.y + 0.15 * (p1.y - my))
    w2 = Point(p1.x + 0.55 * (mx - p1.x), p1.y + 0.55 * (my - p1.y))
    return [InteriorCell(w1, 1), InteriorCell(w2, 2)]


def build_auto_spec(
    domain: DomainInput,
    p1: Point | None = None,
    n_flags: int = 2,
    tau: list[int] | None = None,
    use_collinear: bool | None = None,
) -> ProblemSpec:
    p1 = p1 or find_interior_point(domain.vertices, domain.holes)
    tau_outer = tau or auto_story_tau(domain.vertices, p1, n_flags)
    boundary = polygon_boundary(domain.vertices, tau_outer)
    for hole in domain.holes:
        ht = [1] * len(hole)
        boundary.extend(polygon_boundary(hole, ht))
    pin = auto_pin_route_a(domain.vertices, tau_outer)
    interior = auto_interior_witnesses(p1, pin, n_flags)
    if use_collinear is None:
        xmin, xmax, ymin, ymax = bbox(domain.vertices)
        use_collinear = (xmax - xmin) > 1.2 * max(ymax - ymin, 1e-6)
    collinear = p1.y if use_collinear else None
    return ProblemSpec(
        boundary=boundary,
        p1=p1,
        n_flags=n_flags,
        pins=[pin],
        interior=interior,
        collinear_y=collinear,
        boundary_samples_per_edge=10,
    )


# ---------------------------------------------------------------------------
# 解族扫描：从参数化闭式中找 Viol PASS 的实解
# ---------------------------------------------------------------------------


def _sympy_names_for_flags(n: int) -> list[str]:
    names = []
    for i in range(2, n + 1):
        names.extend([f"p{i}x", f"p{i}y", f"r{i}"])
    return names


def instantiate_parametric(
    param: dict[str, str], subs: dict[str, float]
) -> dict[str, float] | None:
    out: dict[str, float] = {}
    for key, text in param.items():
        expr = sp.sympify(text)
        free = {str(s) for s in expr.free_symbols}
        for s in free:
            if s not in subs:
                return None
        val_c = complex(sp.N(expr.subs({sp.Symbol(k): v for k, v in subs.items()})))
        if abs(val_c.imag) > 1e-7:
            return None
        out[key] = val_c.real
    for i in range(2, 10):
        if f"r{i}" in out and out[f"r{i}"] <= 0:
            return None
    return out


def free_symbols_in_param(param: dict[str, str]) -> set[str]:
    free: set[str] = set()
    for text in param.values():
        free |= {str(s) for s in sp.sympify(text).free_symbols}
    return free


def score_solution(spec: ProblemSpec, sol: dict[str, float]) -> float:
    pos, rad = positions_radii_from_solution(spec, sol)
    viol_ok, viol_logs = verify_viol(spec, pos, rad, samples_per_edge=12)
    int_ok, _ = verify_interior(spec, pos, rad)
    if viol_ok and int_ok:
        return 0.0
    return 10.0 + len(viol_logs) + (0 if int_ok else 5.0)


def aesthetic_score(spec: ProblemSpec, sol: dict[str, float]) -> float:
    """Viol 同为 0 时，优先选共线、非贴边、半径合理的代表点。"""
    pen = 0.0
    px = sol.get("p2x", 0.0)
    py = sol.get("p2y", 0.0)
    r2 = sol.get("r2", 0.0)
    if spec.collinear_y is not None:
        pen += 80.0 * abs(py - spec.collinear_y)
    pen += max(0.0, 0.18 - py) * 120.0
    pen += max(0.0, px - spec.p1.x - 1.8) * 10.0
    if px < spec.p1.x + 0.04:
        pen += 25.0
    if r2 < 0.2:
        pen += 30.0
    if r2 > 2.5:
        pen += 15.0
    return pen


def combined_score(spec: ProblemSpec, sol: dict[str, float]) -> float:
    return score_solution(spec, sol) + aesthetic_score(spec, sol)


def resolve_solution(
    spec: ProblemSpec, report: SolveReport
) -> tuple[dict[str, float] | None, str]:
    """从消元报告挑选或扫描出可绘制的实解。"""
    candidates: list[dict[str, float]] = []
    notes: list[str] = []

    for sol in report.solutions:
        if score_solution(spec, sol) < 1e-9:
            return sol, "finite 实解"
        candidates.append(sol)

    if report.parametric:
        for fam in report.parametric:
            free = sorted(free_symbols_in_param(fam))
            if not free:
                continue
            if len(free) == 1:
                name = free[0]
                lo, hi = 0.05, 4.0
                if name.startswith("r"):
                    lo, hi = 0.08, 3.5
                best: dict[str, float] | None = None
                best_score = 1e9
                for val in np.linspace(lo, hi, 120):
                    inst = instantiate_parametric(fam, {name: float(val)})
                    if inst is None:
                        continue
                    sc = score_solution(spec, inst)
                    if sc < best_score:
                        best_score = sc
                        best = inst
                    if sc < 1e-9:
                        return inst, f"参数族扫描 ({name}={val:.4f})"
                if best is not None and best_score < 5.0:
                    candidates.append(best)
                    notes.append(f"参数族近似 ({name})")
            elif len(free) == 2:
                n1, n2 = free[0], free[1]
                best = None
                best_score = 1e9
                y_center = spec.collinear_y if spec.collinear_y is not None else 0.35
                if "p2y" in free:
                    y_grid = np.linspace(max(0.05, y_center - 0.2), y_center + 0.2, 31)
                else:
                    y_grid = np.linspace(0.08, 2.5, 25)
                x_grid = np.linspace(0.1, 2.5, 25)
                r_grid = np.linspace(0.2, 2.5, 25)
                grids = {
                    "p2x": x_grid,
                    "p2y": y_grid,
                    "r2": r_grid,
                }
                g1 = grids.get(n1, np.linspace(0.1, 2.5, 25))
                g2 = grids.get(n2, np.linspace(0.08, 2.5, 25))
                for v1 in g1:
                    for v2 in g2:
                        inst = instantiate_parametric(fam, {n1: float(v1), n2: float(v2)})
                        if inst is None:
                            continue
                        sc = combined_score(spec, inst)
                        if sc < best_score:
                            best_score = sc
                            best = inst
                if best is not None and score_solution(spec, best) < 8.0:
                    candidates.append(best)
                    notes.append(f"二维参数族 ({n1},{n2})")

    if candidates:
        best = min(candidates, key=lambda s: combined_score(spec, s))
        note = notes[0] if notes else "最优候选（Viol 最小）"
        return best, note
    return None, "无可行实解"


# ---------------------------------------------------------------------------
# 可视化
# ---------------------------------------------------------------------------

FLAG_COLORS = {
    1: "#4C72B0",
    2: "#DD8452",
    3: "#55A868",
    4: "#C44E52",
    5: "#8172B3",
    6: "#937860",
}


def compute_region_grid(
    outer: Sequence[tuple[float, float]],
    holes: Sequence[Sequence[tuple[float, float]]],
    positions: list[tuple[float, float]],
    radii: list[float],
    res: int = 320,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xmin, xmax, ymin, ymax = bbox(outer, 0.02)
    xs = np.linspace(xmin, xmax, res)
    ys = np.linspace(ymin, ymax, res)
    X, Y = np.meshgrid(xs, ys)
    labels = np.zeros_like(X, dtype=np.int32)
    mask = np.zeros_like(X, dtype=bool)
    for iy in range(res):
        for ix in range(res):
            x, y = float(X[iy, ix]), float(Y[iy, ix])
            if point_in_domain(x, y, outer, holes):
                mask[iy, ix] = True
                labels[iy, ix] = assign_min_flag(x, y, positions, radii)
    return X, Y, labels, mask


def plot_tyson(
    domain: DomainInput,
    spec: ProblemSpec,
    sol: dict[str, float],
    out_path: Path,
    dpi: int = 150,
    grid_res: int = 360,
) -> None:
    pos, rad = positions_radii_from_solution(spec, sol)
    n = spec.n_flags
    xmin, xmax, ymin, ymax = bbox(domain.vertices, 0.08)
    fig, ax = plt.subplots(figsize=(9, 8))

    X, Y, labels, mask = compute_region_grid(
        domain.vertices, domain.holes, pos, rad, res=grid_res
    )

    # 势力区填色
    for flag in range(1, n + 1):
        c = FLAG_COLORS.get(flag, "#888888")
        region = np.ma.masked_where((labels != flag) | ~mask, labels)
        ax.contourf(X, Y, region, levels=[flag - 0.5, flag + 0.5], colors=[c], alpha=0.42)

    # 域边界
    ox, oy = zip(*domain.vertices)
    ax.plot(list(ox) + [ox[0]], list(oy) + [oy[0]], "k-", lw=2.2, label="boundary")
    for hole in domain.holes:
        hx, hy = zip(*hole)
        ax.plot(list(hx) + [hx[0]], list(hy) + [hy[0]], "k--", lw=1.8)

    # 圆盘 B(pi, ri)
    for i in range(n):
        cx, cy = pos[i]
        r = rad[i]
        circ = Circle((cx, cy), r, fill=False, linestyle="--", linewidth=1.2,
                      edgecolor=FLAG_COLORS.get(i + 1, "#333"), alpha=0.85)
        ax.add_patch(circ)

    # 站点
    for i in range(n):
        cx, cy = pos[i]
        ax.scatter([cx], [cy], s=120, c=FLAG_COLORS.get(i + 1, "#333"), edgecolors="k", zorder=5)
        ax.annotate(
            f"$p_{i+1}$", (cx, cy), textcoords="offset points", xytext=(8, 8),
            fontsize=12, fontweight="bold",
        )

    # 钉点
    for pin in spec.pins:
        ax.scatter(
            [pin.point.x], [pin.point.y], marker="x", s=80, c="red", linewidths=2, zorder=6
        )
        ax.annotate(
            f"Π∈Γ_{pin.i}{pin.j}",
            (pin.point.x, pin.point.y),
            textcoords="offset points",
            xytext=(6, -14),
            fontsize=9,
            color="red",
        )

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    title = "Mul-Tyson  " + ", ".join(
        f"p{i+1}=({pos[i][0]:.3f},{pos[i][1]:.3f}), r{i+1}={rad[i]:.3f}"
        for i in range(n)
    )
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def domain_from_spec_json(path: Path) -> DomainInput:
    import json

    data = json.loads(path.read_text(encoding="utf-8"))
    verts = [tuple(v) for v in data["vertices"]]
    holes = [list(map(tuple, h)) for h in data.get("holes", [])]
    return DomainInput(vertices=verts, holes=holes)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="任意域乘性泰森：求解 + 可视化")
    parser.add_argument("--vertices", help='外边界顶点 "x,y x,y ..."')
    parser.add_argument("--holes", action="append", help='孔洞顶点（可多次）')
    parser.add_argument("--spec", help="JSON 规格（含 vertices）")
    parser.add_argument("--p1", help="起点 p1；省略则自动找内点")
    parser.add_argument("--n", type=int, default=2, help="旗数 N")
    parser.add_argument("--tau", help="边归属 edge:tau，如 0:1 1:2")
    parser.add_argument("--pin", action="append", help="手动硬钉 x,y,i,j")
    parser.add_argument("--no-auto-story", action="store_true", help="不自动设 τ/钉点")
    parser.add_argument("--out", default="tyson_plot.png", help="输出图片路径")
    parser.add_argument("--grid", type=int, default=360, help="栅格分辨率")
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--report", action="store_true", help="打印消元报告")
    args = parser.parse_args(argv)

    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    # 域
    if args.spec:
        domain = domain_from_spec_json(Path(args.spec))
        spec = load_spec(Path(args.spec))
    elif args.vertices:
        verts = parse_vertices(args.vertices)
        holes = [parse_vertices(h) for h in (args.holes or [])]
        domain = DomainInput(vertices=verts, holes=holes)
        p1 = None
        if args.p1:
            px, py = map(float, args.p1.split(","))
            p1 = Point(px, py)
        tau = parse_tau(args.tau or "", len(verts)) if args.tau else None
        if args.no_auto_story and args.pin:
            pins = []
            for pin_s in args.pin:
                xs, ys, i_s, j_s = pin_s.split(",")
                pins.append(Pin(Point(float(xs), float(ys)), int(i_s), int(j_s)))
            spec = ProblemSpec(
                boundary=polygon_boundary(verts, tau or [1] * len(verts)),
                p1=p1 or find_interior_point(verts, holes),
                n_flags=args.n,
                pins=pins,
            )
        else:
            spec = build_auto_spec(domain, p1=p1, n_flags=args.n, tau=tau)
            if args.pin:
                for pin_s in args.pin:
                    xs, ys, i_s, j_s = pin_s.split(",")
                    spec.pins.append(Pin(Point(float(xs), float(ys)), int(i_s), int(j_s)))
    else:
        parser.error("需要 --vertices 或 --spec")

    print(f"p1 = ({spec.p1.x:.4f}, {spec.p1.y:.4f})")
    print(f"自动故事: |Π|={len(spec.pins)}, τ 边数={len(spec.boundary)}")

    report = solve_spec(spec)
    if args.report:
        print(format_report(spec, report))

    sol, note = resolve_solution(spec, report)
    if sol is None:
        print("未找到可可视化的实解。尝试 --report 查看消元结果，或手动 --pin / --tau。")
        return 1

    print(f"选用解: {note}")
    for i in range(2, spec.n_flags + 1):
        print(f"  p{i}=({sol[f'p{i}x']:.6f}, {sol[f'p{i}y']:.6f}), r{i}={sol[f'r{i}']:.6f}")

    out = Path(args.out)
    plot_tyson(domain, spec, sol, out, dpi=args.dpi, grid_res=args.grid)
    print(f"已保存: {out.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
