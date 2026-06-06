#!/usr/bin/env python3
"""验证路线 C 守卫参数：矩形 §21.1 与 L 形 §21.2（非共点四旗）。

非证明链：主解析链见 symmetric_sector_theorems.md §21；本脚本仅复算。

用法:
  python verify_guard_params.py
  python verify_guard_params.py --shape rect
  python verify_guard_params.py --shape L --xi 0.55 --eps 0.1 --r 0.9
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass


@dataclass
class Config:
    name: str
    W: float  # Rx for L
    H: float  # w for L (arm height)
    a: float
    b: float
    xi: float
    eps: float
    r: float
    x_leak_min: float  # only audit leakage segment on horizontal edges


def chi12(x: float, y: float, a: float, b: float, mid: float, lam: float) -> float:
    d1 = math.hypot(x - a, y - mid)
    d2 = math.hypot(x - b, y - mid) / lam
    return d2 - d1


def region(
    x: float, y: float, a: float, b: float, mid: float, lam: float, p3, r3, p4, r4
) -> int:
    d1 = math.hypot(x - a, y - mid)
    d2 = math.hypot(x - b, y - mid) / lam
    d3 = math.hypot(x - p3[0], y - p3[1]) / r3
    d4 = math.hypot(x - p4[0], y - p4[1]) / r4
    dists = [(1, d1), (2, d2), (3, d3), (4, d4)]
    winner = min(dists, key=lambda t: t[1])[0]
    return winner


def verify(cfg: Config, n_sample: int = 400) -> bool:
    mid = cfg.H / 2
    lam = (cfg.W - cfg.b) / (cfg.W - cfg.a)
    dM = cfg.W - cfg.a
    M = (cfg.W, mid)

    p3 = (cfg.xi, cfg.H - cfg.eps)
    p4 = (cfg.xi, cfg.eps)

    ok = True
    lines: list[str] = []

    def log(msg: str) -> None:
        lines.append(msg)

    log(f"=== {cfg.name} ===")
    log(f"p1=({cfg.a},{mid}), p2=({cfg.b},{mid}), lambda={lam:.6f}")
    log(f"p3={p3}, p4={p4}, r3=r4={cfg.r}")

    for label, p in [("M vs p3", p3), ("M vs p4", p4)]:
        d1 = math.hypot(M[0] - cfg.a, M[1] - mid)
        d_guard = math.hypot(M[0] - p[0], M[1] - p[1]) / cfg.r
        safe = d1 <= d_guard + 1e-9
        log(f"  {label}: d(M,p1)={d1:.4f}, d(M,p)/r={d_guard:.4f}  M_safe={safe}")
        ok = ok and safe

    xs_top = [
        cfg.x_leak_min + (cfg.W - cfg.x_leak_min) * i / (n_sample - 1)
        for i in range(n_sample)
    ]
    xs_bot = list(xs_top)

    v2_count = 0
    total = 0
    for y, pname in [(cfg.H, "top"), (0.0, "bottom")]:
        for x in xs_top:
            d1 = math.hypot(x - cfg.a, y - mid)
            d2 = math.hypot(x - cfg.b, y - mid) / lam
            d3 = math.hypot(x - p3[0], y - p3[1]) / cfg.r
            d4 = math.hypot(x - p4[0], y - p4[1]) / cfg.r
            d_min_other = min(d1, d3, d4)
            if d2 < d_min_other - 1e-7:
                v2_count += 1
            total += 1

    log(f"  horizontal edges: V2 points = {v2_count}/{total}")
    ok = ok and v2_count == 0

    log(f"RESULT: {'PASS' if ok else 'FAIL'}")
    print("\n".join(lines))
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify guard route-C parameters")
    parser.add_argument("--shape", choices=["rect", "L", "both"], default="both")
    parser.add_argument("--xi", type=float)
    parser.add_argument("--eps", type=float)
    parser.add_argument("--r", type=float)
    args = parser.parse_args()

    presets = {
        "rect": Config(
            name="Rectangle (Sec 21.1)",
            W=2.0,
            H=1.0,
            a=0.5,
            b=0.7,
            xi=0.9,
            eps=0.15,
            r=0.7,
            x_leak_min=0.82,
        ),
        "L": Config(
            name="L-shape horizontal arm (Sec 21.2)",
            W=2.0,
            H=0.5,
            a=0.4,
            b=0.6,
            xi=0.55,
            eps=0.1,
            r=0.9,
            x_leak_min=0.55,
        ),
    }

    shapes = ["rect", "L"] if args.shape == "both" else [args.shape]
    all_ok = True
    for key in shapes:
        cfg = presets[key]
        if args.xi is not None:
            cfg.xi = args.xi
        if args.eps is not None:
            cfg.eps = args.eps
        if args.r is not None:
            cfg.r = args.r
        all_ok = verify(cfg) and all_ok

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
