#!/usr/bin/env python3
"""算法 22.1 / 定理 M★ 审计（§23.5–23.6 工具层，非证明链）。

在特设 Σ + 强规范 B* + 路线 A 下，对文档 §18.7 登记算例执行：
  步骤 2：m_eff=0 ⇒ n_min=2
  步骤 3：m_eff≥1 ⇒ 秩3 (p3,r3)=(p2,λ) ⇒ n_min=3

用法:
  python nmin_decision_audit.py
  python nmin_decision_audit.py --verbose
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass


@dataclass
class Instance:
    id: str
    m_eff: int
    n_min_rank3: int
    n_min_geo: int | None  # None = 不适用
    source: str
    note: str = ""


# §18.7 / §23.4 登记实例（与形状名无关的审计结论）
INSTANCES: list[Instance] = [
    Instance("solid_sector", 0, 2, None, "§8 M.8", "OA,OB 无渗漏"),
    Instance("concentric_hole_sector", 0, 2, None, "§9 M.9", "孔弧全 χ≥0"),
    Instance("eccentric_hole_arc", 1, 3, None, "§11.6 M.11", "秩3 定理 19.4"),
    Instance("rect_top_only", 1, 3, None, "§19 M.19", "单顶边 m_eff=1"),
    Instance("rect_top_bottom", 2, 3, 4, "§21 M.21", "秩3=3；geo 双钉=4"),
    Instance("L_boundary_Bstar", 2, 3, 4, "§21 M.15", "秩3=3；geo 定理 21.2=4"),
    Instance("L_vertical_arm", 0, 2, None, "§15", "竖臂验收 m_eff=0"),
]


def rank3_no_v2_strict_global() -> tuple[bool, str]:
    """引理 M★.0 数值烟测：p3=p2,r3=λ 时 φ2≡φ3，无 V2 严格内点。"""
    a, b, lam = 0.5, 0.7, 13 / 15
    p2 = (b, 0.5)
    p3 = p2
    r3 = lam
    mid = 0.5
    strict_v2 = 0
    for i in range(50):
        for j in range(50):
            x = 0.1 + 1.8 * i / 49
            y = 0.05 + 0.9 * j / 49
            d1 = math.hypot(x - a, y - mid)
            d2 = math.hypot(x - b, y - mid) / lam
            d3 = math.hypot(x - p3[0], y - p3[1]) / r3
            if d2 < d1 - 1e-9 and d2 < d3 - 1e-9:
                strict_v2 += 1
    ok = strict_v2 == 0
    msg = f"rank3 global: V2 strict interior count = {strict_v2} (expect 0)"
    return ok, msg


def audit_algorithm_22_1() -> list[tuple[Instance, str, bool]]:
    """模拟算法 22.1 秩3 分支（步骤 2 / 3）。"""
    rows: list[tuple[Instance, str, bool]] = []
    for inst in INSTANCES:
        if inst.m_eff == 0:
            step = 2
            n_out = 2
            branch = "M.I"
        else:
            step = 3
            n_out = 3
            branch = "M.II 秩3"
        ok = n_out == inst.n_min_rank3
        detail = f"step{step} {branch} → n_min={n_out}"
        if inst.n_min_geo is not None:
            detail += f" | n_min^geo={inst.n_min_geo}（几何分离，非秩3）"
        rows.append((inst, detail, ok))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="算法 22.1 / 定理 M★ 审计")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    print("=" * 60)
    print("定理 M★ 审计：Σ + B* + 路线 A，秩3 路线 n_min ∈ {2, 3}")
    print("=" * 60)

    ok_global, msg_global = rank3_no_v2_strict_global()
    print(f"[{'PASS' if ok_global else 'FAIL'}] {msg_global}")

    all_ok = ok_global
    print()
    print(f"{'实例':<28} {'m_eff':>5} {'期望':>4}  算法 22.1 输出")
    print("-" * 60)
    for inst, detail, ok in audit_algorithm_22_1():
        all_ok = all_ok and ok
        mark = "PASS" if ok else "FAIL"
        print(
            f"[{mark}] {inst.id:<24} {inst.m_eff:>5} {inst.n_min_rank3:>4}  {detail}"
        )
        if args.verbose:
            print(f"       来源 {inst.source}: {inst.note}")

    print("-" * 60)
    print("未闭合分支（仅 n_min^geo 几何分离路线，§23.8）：")
    print("  - M.IV-E：非共点 14.4 守卫存在性（m_eff≥2 时 n_min^geo=2+m_eff）")
    print("  - M.V+：轴上分离等障碍族对 geo 路线是否完备")
    print("  - 孔弧 p3!=p2 其它几何分离（命题 19.6 仍开放）")
    print()
    print(f"OVERALL: {'PASS' if all_ok else 'FAIL'}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
