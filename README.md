# 乘性泰森地皮划分

## 问题

给定平面区域 \(\Omega\)、第一面旗位置 \(p_1\)、边界归属 \(\tau\)（每条边归哪面旗），反求其余旗位与半径，并判定解空间：

- 无解
- 有限多解（全列）
- 无穷多解（参数化闭式）

同时按定理 M★ 自动给出 \(n_{\min}\)（\(m_{\mathrm{eff}}=0\Rightarrow 2\)，\(m_{\mathrm{eff}}\ge 1\Rightarrow 3\) 秩3路线）。

## 本仓库内容

数值复现与审计脚本，对应扇形 Σ9 分类、宽题 W 封口（§W.12）、路线 C 守卫参数等。**非证明链**；主文稿另行整理。

## 依赖

```bash
pip install -r requirements.txt
```

## 主要脚本

| 脚本 | 用途 |
|------|------|
| `solve_land.py` | 原问题封口：\(\Omega,p_1,\tau\) → \(n_{\min}\) + 解空间报告 |
| `mul_tyson_solve.py` | Gröbner 消元核心 |
| `mul_tyson_viz.py` | 求解 + 势力区可视化 |
| `sector_sigma9_verify.py` | 扇形 Σ9 符号/数值验证 |
| `nmin_decision_audit.py` | 定理 M★ 算例审计 |
| `verify_guard_params.py` | 矩形/L 形路线 C 守卫参数复算 |
| `decidability_convex_pin.py` | 凸域单钉自由度实验 |
| `compare_boundary_story.py` | 与 CGAL/容量泰森对照 |
| `plot_sector_figures.py` / `plot_l_shape_figures.py` | 标准算例示意图生成（本地 `-o` 输出） |

## 示例

```bash
python solve_land.py --vertices "0,0 2,0 2,1 0,1" --p1 "0.4,0.5" --tau "0:1 1:2 2:1 3:1"
python sector_sigma9_verify.py
python nmin_decision_audit.py
python verify_guard_params.py
```
