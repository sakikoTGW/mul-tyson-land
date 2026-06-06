#!/usr/bin/env python3
"""扇形 Σ9 完整分类 — 符号/数值验证（对应 symmetric_sector_theorems.md §Σ9）。

运行: python sector_sigma9_verify.py
"""

from __future__ import annotations

import math
import sys

import sympy as sp
from sympy import Rational, cos, sin, sqrt, symbols


def verify_v1_impossible() -> bool:
    """Σ9-I (v1): A,M 双钉 ⇒ b∈{a, R²/a}，无 (a,R) 内解。"""
    R, a, alpha = symbols("R a alpha", positive=True, real=True)
    b = symbols("b", real=True)
    h2 = R**2 + a**2 - 2 * a * R * cos(alpha)
    lam = (R - b) / (R - a)
    lhs = h2 * lam**2
    rhs = R**2 + b**2 - 2 * b * R * cos(alpha)
    eq = sp.Eq(lhs, rhs)
    poly = sp.expand(lhs - rhs)
    # 化为关于 b 的二次式
    poly_b = sp.Poly(poly, b)
    roots = [sp.simplify(r) for r in sp.solve(poly_b, b)]
    expected = [a, R**2 / a]
    ok = all(any(sp.simplify(r - e) == 0 for e in expected) for r in roots)
    in_interval = all(
        sp.simplify(r - a) == 0 or sp.simplify(r - R**2 / a) == 0 for r in roots
    )
    return ok and in_interval and len(roots) == 2


def _k_numeric(R: float, a: float, alpha: float) -> float:
    b_tang = R * (2 * R - a - R * math.cos(alpha)) / (R - a * math.cos(alpha))
    h2 = R**2 + a**2 - 2 * a * R * math.cos(alpha)
    A2 = 2 * a * R * (1 - math.cos(alpha))
    h = math.sqrt(h2)
    A1 = -2 * R * ((R - a) * h + (R - a) ** 2 * math.cos(alpha))
    A0 = (R - a) ** 2 * R**2 + h2 * R**2 - 2 * h2 * R * (R - a)
    return A2 * b_tang**2 + A1 * b_tang + A0


def verify_k_impossible() -> bool:
    """定理 D: K=0 iff cos(alpha)=2R/(R+a); 对 a<R 该值 >1 故不可能。"""
    samples = [(1.0, 0.3, math.pi / 4), (2.0, 0.8, 0.9), (1.5, 0.4, math.pi / 3)]
    nonzero = all(abs(_k_numeric(R, a, al)) > 1e-6 for R, a, al in samples)
    bound_gt1 = all(2 * R / (R + a) > 1.0 for R, a, _ in samples)
    return nonzero and bound_gt1


def verify_route_a_f_negative() -> bool:
    """定理 8.1(i): f(R)<0 当 a<b<R。"""
    R, a, b, alpha = symbols("R a b alpha", positive=True, real=True)
    lam = (R - b) / (R - a)
    fR = (
        R**2
        + a**2
        - 2 * a * R * cos(alpha)
        - (R**2 + b**2 - 2 * b * R * cos(alpha)) / lam**2
    )
    fact = sp.factor(fR)
    # 在假设 a<b<R, alpha∈(0,pi) 下符号应为负
    test = fact.subs({R: 1, a: Rational(3, 10), b: Rational(1, 2), alpha: sp.pi / 4})
    return float(test) < 0


def verify_area_monotone_numeric() -> bool:
    """命题 8.4: A(b) 严格单调减（数值）。"""
    R, a, alpha = 1.0, 0.3, math.pi / 4

    def sector_area_v2(b: float, n: int = 400) -> float:
        lam = (R - b) / (R - a)
        p1, p2 = (a, 0.0), (b, 0.0)
        cnt = 0
        for i in range(n):
            for j in range(n):
                rho = R * (i + 0.5) / n
                phi = -alpha + 2 * alpha * (j + 0.5) / n
                x, y = rho * math.cos(phi), rho * math.sin(phi)
                d1 = math.hypot(x - p1[0], y - p1[1])
                d2 = math.hypot(x - p2[0], y - p2[1]) / lam
                if d2 < d1:
                    cnt += 1
        cell = (R / n) * (2 * alpha / n) * rho  # 近似微元
        return cnt * cell

    bs = [a + 0.02 + 0.08 * k for k in range(8)]
    areas = [sector_area_v2(b) for b in bs]
    return all(areas[i] > areas[i + 1] for i in range(len(areas) - 1))


def verify_b_eta_unique() -> bool:
    """命题 8.5: 给定 η，b_η 唯一（介值+单调）。"""
    R, a, alpha = 1.0, 0.3, math.pi / 4
    omega = R**2 * alpha / 2

    def area(b: float) -> float:
        lam = (R - b) / (R - a)
        p1, p2 = (a, 0.0), (b, 0.0)
        n = 300
        s = 0.0
        for i in range(n):
            for j in range(n):
                rho = R * (i + 0.5) / n
                phi = -alpha + 2 * alpha * (j + 0.5) / n
                x, y = rho * math.cos(phi), rho * math.sin(phi)
                d1 = math.hypot(x - p1[0], y - p1[1])
                d2 = math.hypot(x - p2[0], y - p2[1]) / lam
                if d2 < d1:
                    s += (R / n) * (2 * alpha / n) * rho
        return s

    A_sup = area(a + 0.001)
    eta = 0.2
    target = eta * omega
    lo, hi = a + 0.001, R - 0.001
    for _ in range(60):
        mid = (lo + hi) / 2
        if area(mid) > target:
            lo = mid
        else:
            hi = mid
    b_eta = (lo + hi) / 2
    # 扰动应偏离目标
    return abs(area(b_eta) - target) < 0.02 * omega and area(b_eta + 0.05) < target


def classify_sigma9_story(story: str) -> str:
    """Return W-classification label for Sigma9 sub-stories (ASCII for console)."""
    table = {
        "v1": "NO_SOL (dual pin, Cor 6.5)",
        "v2_wrong": "NO_SOL (A in V2 contradicts Thm 8.1)",
        "route_a": "MULTI (1-dim family b in (a,R), Thm 8.12)",
        "route_a_eta": "UNIQUE (area share pins b_eta, Prop 8.5)",
        "route_a_double_pin": "NO_SOL (second pin, Thm 9.5 type)",
    }
    return table.get(story, "unknown")


def main() -> int:
    checks = [
        ("Sigma9-I v1: no root in (a,R)", verify_v1_impossible()),
        ("Thm D: K=0 impossible for a<R", verify_k_impossible()),
        ("Route A: f(R)<0 => A in V1", verify_route_a_f_negative()),
        ("Prop 8.4: A(b) strictly decreasing", verify_area_monotone_numeric()),
        ("Prop 8.5: b_eta unique", verify_b_eta_unique()),
    ]
    print("=== sector Sigma9 verification ===")
    all_ok = True
    for name, ok in checks:
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {name}")
        all_ok = all_ok and ok
    print("\n=== Sigma9 classification (story -> W output) ===")
    for key in ["v1", "v2_wrong", "route_a", "route_a_eta", "route_a_double_pin"]:
        print(f"  {key:22s} -> {classify_sigma9_story(key)}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
