# 乘性泰森地皮划分

## 1　引言

设 $`\Omega\subset\mathbb R^2`$ 为有界闭域，$`\partial\Omega`$ 由分段 $C^1$ 折线围成。$`\Omega`$ 的边界可为任意多边形（含海星形等非凸外形）；亦允许孔洞，例如环域

$$
\partial\Omega=\Gamma_{\mathrm{out}}\sqcup\bigsqcup_h\Gamma_{\mathrm{in}}^{(h)}.
$$

外圈为 $`\Gamma_{\mathrm{out}}`$，内圈为孔边界 $`\Gamma_{\mathrm{in}}^{(h)}`$。给定 $`p_1\in\Omega`$，第一面旗半径归一 $`r_1=1`$。对其余旗 $`(p_i,r_i)`$（$`r_i>0`$），定义

$$
\varphi_i(x):=\frac{\|x-p_i\|}{r_i},
$$

势力区定义为

$$
V_i=\{x\in\Omega:\varphi_i(x)\le\varphi_j(x),\,\forall j\neq i\}.
$$

旗 $i$ 与 $j$ 的分界线 $`\{\varphi_i=\varphi_j\}`$ 为阿波罗尼奥斯圆；$`r_i=r_j`$ 时退化为直线。边界标签 $`\tau`$ 指定 $`\partial\Omega`$ 各边（除分界点外）的目标归属。要求

$$
\bigcup_{i=1}^N V_i=\Omega.
$$

在边界 $`\partial\Omega`$ 上（分界点除外），还要求归属一致：

$$
\sigma_{\mathcal S}(s)=\tau(s),\qquad s\in\partial\Omega,
$$

其中

$$
\sigma_{\mathcal S}(x)=\arg\min_{i}\varphi_i(x).
$$

**反问题**　已知 $`\Omega,p_1,\tau`$，求 $`(p_2,\ldots,p_N,r_2,\ldots,r_N)`$，并判定 $`\mathfrak{Sol}`$ 为无解、唯一、有限或无穷；无穷时给出参数化闭式。

在路线 A（(1.2)(1.3)）与强规范 $`\mathcal B^{\ast}`$（渗漏边集 $`\mathcal E=\varnothing`$，见 §4）下，记

$$
m_{\mathrm{eff}}=\lvert\mathcal E\rvert.
$$

**定理 1.1**

$$
m_{\mathrm{eff}}=0\Rightarrow n_{\min}=2;\qquad m_{\mathrm{eff}}\ge1\Rightarrow n_{\min}=3,\quad (p_3,r_3)=(p_2,\lambda).
$$

**定理 1.2**

$$
\Omega=[0,W]\times[0,H],\quad p_1=(a,H/2),\quad M=(W,H/2)\in\Gamma_{12}
\ \Longrightarrow\ b=W-(W-a)r_2\quad (r_2>0),
$$

为解集的一条 1 维族。

**定理 1.3**　扇形 $`\Omega(R,\alpha)`$、路线 A 的 1 维族上，若

$$
\mathrm{Area}(V_2)=\eta\lvert\Omega\rvert,\qquad \eta\in(0,1),
$$

则存在唯一 $`b_\eta\in(a,R)`$。

---

## 2　几何对象与划分

### §1.1 地皮

地皮是一块平面**闭域** $`\Omega\subset\mathbb R^2`$：

- 外边界 $`\Gamma_{\mathrm{out}}`$ 可由折线围成任意多边形（矩形、L 形、扇形、海星形等均为**代入同一框架的算例**，不另立形状定理）；
- 允许**孔洞**。边界分解为

  $$
  \partial\Omega=\Gamma_{\mathrm{out}}\sqcup\bigsqcup_h\Gamma_{\mathrm{in}}^{(h)}.
  $$

  孔边界与外边界的归属均由 $`\tau`$ 指定；
- 假设边界分段 $C^1$，除有限个角点外法向可定义。

### §1.2 旗与性价比

在平地上选定一点插**第一面旗** $`p_1\in\Omega`$，管辖半径取 $`r_1=1`$（归一化，不失一般性）。

尚需插入 $N-1$ 面旗。第 $i$ 面旗有位置 $`p_i\in\mathbb R^2`$ 与管辖半径 $`r_i>0`$。配置记为

$$
\mathcal S=\{p_1,\ldots,p_N\}\cup\{r_1,\ldots,r_N\},\qquad r_1=1.
$$

对 $x\in\mathbb R^2$，到第 $i$ 面旗的**性价比**（乘性泰森权）定义为

$$
\boxed{\varphi_i(x):=\frac{\|x-p_i\|}{r_i}} \qquad (1.1)
$$

$x$ 归 $`\varphi_i`$ **最小**的那面旗管辖。记势力区

$$
V_i:=\{x\in\Omega:\ \varphi_i(x)\le\varphi_j(x)\ \forall j\neq i\}
$$

（分界线上归属可任取一侧，零测度不影响面积与审计）。

**支配场**：

$$
\sigma_{\mathcal S}(x):=\arg\min_{1\le i\le N}\varphi_i(x).
$$

### §1.3 势力分界线：阿波罗尼奥斯圆

两旗 $`p_i,p_j`$ 性价比相等时

$$
\frac{\|x-p_i\|}{r_i}=\frac{\|x-p_j\|}{r_j}
\quad\Longleftrightarrow\quad
\|x-p_i\|^2 r_j^2=\|x-p_j\|^2 r_i^2.
$$

这是平面上的圆（阿波罗尼奥斯圆）；仅当 $`r_i=r_j`$ 时退化为垂直平分直线。Laguerre 幂图 $`\pi_i=\|x-p_i\|^2-r_i^2`$ 的胞间分界为直线；本文采用乘性泰森 $`\varphi_i=d/r_i`$。

### §1.4 边界故事 $`\tau`$ 与可行性

**边界故事** $`\tau`$：对 $`\partial\Omega`$ 的每条边（或子弧）$e$，指定 $e$ 上除分界点外的目标归属 $`\tau(e)\in\{1,\ldots,N\}`$。

**整包约束** $`\mathcal G`$：

| 编号 | 约束 | 形式 |
|------|------|------|
| G1 覆盖 | 整块地分完 | $`V_1\cup\cdots\cup V_N=\Omega`$ |
| G2 缝相容 | 拼接片在缝上归属一致 | 分片 $`\Omega`$ 时边界标签连续 |
| G3 不过定 | 钉点与自由未知数相容 | 见 §6 |
| $`\mathcal B`$ | 边界故事成立 | $`\sigma_{\mathcal S}(s)=\tau(s)`$ 对 $`s\in\partial\Omega`$（分界点除外） |

**违约集**：

$$
\mathrm{Viol}(\mathcal S;\tau):=\{s\in\partial\Omega:\ \sigma_{\mathcal S}(s)\neq\tau(s)\}.
$$

**引理 2.1**　$`\mathfrak{Sol}`$ 为半代数集。若 $E$ 个钉方程在一般位置独立，则 $`\dim\mathfrak{Sol}\ge\max(0,\mathrm{dof}-E)`$。

证明. 每个独立多项式等式至多使解集降维 1。若 $`E<\mathrm{dof}`$，维数严格为正。证毕.

**引理 2.2**　$`r_1\neq r_2`$ 时 $`\{\varphi_1=\varphi_2\}`$ 为圆；$`r_1=r_2`$ 时为直线。

证明. $`\|x-p_1\|/r_1=\|x-p_2\|/r_2`$ 平方整理即 $`\|x-p_1\|^2 r_2^2=\|x-p_2\|^2 r_1^2`$，为标准阿波罗尼奥斯圆方程。证毕.

**引理 2.3**　$N=2$ 时 $`V_1\cup V_2=\Omega`$。

证明. 对 $`x\in\Omega`$，$`\varphi_1(x)\le\varphi_2(x)`$ 或 $`\varphi_2(x)\le\varphi_1(x)`$。证毕.

---

## 3　反问题

### §2.1 输入与输出

**输入**

1. 地皮 $`\Omega`$（顶点序列 + 可选孔洞）；
2. 第一面旗 $`p_1\in\Omega`$，$`r_1=1`$；
3. 边界故事 $`\tau`$（每条边归哪面旗）。

**过程**

4. 自动判定最少旗数 $`n_{\min}`$（定理 1.1）；
5. 在固定 $N$ 下建立方程组并消元。

**输出**

| 分类 | 数学含义 |
|------|----------|
| **无解** | $`\mathfrak{Sol}=\varnothing`$ 或无可行实解满足 Viol |
| **唯一** | $`\dim\mathfrak{Sol}=0`$ 且 $`\lvert\mathfrak{Sol}\cap\mathfrak D^c\rvert=1`$ |
| **有限多解** | $`\dim\mathfrak{Sol}=0`$ 且 $`\lvert\mathfrak{Sol}\rvert>1`$，全离散 |
| **无穷多解** | $`\dim\mathfrak{Sol}\ge 1`$，给出参数化闭式与自由参 |

其中 $`\mathfrak D`$ 为退化簇（$`p_i=p_j`$、$`r_i=0`$ 等）。

### §2.2 解析解的含义

**定义 3.1（解析解）**　边界匹配代数方程组的实根，或 Gröbner 消元给出的参数化闭式；允许正维解族。

---

## 4　宽题与特设 $\Sigma$

### 4.1　宽题 W

**输入**：有界域 $`\Omega\subset\mathbb R^2`$（任意封闭图形，允许多孔洞 $`\partial\Omega=\Gamma_{\mathrm{out}}\sqcup\bigsqcup_h\Gamma_{\mathrm{in}}^{(h)}`$）；固定 $`p_1\in\Omega`$，$`r_1=1`$；边界故事 $`\mathcal B`$（含弱/强规范）。

**要求**：求 $`(p_2,\ldots,p_N)`$ 与 $`(r_2,\ldots,r_N)>0`$，使势力区 $`V_i`$ 同时满足整包约束 $`\mathcal G`$：

| 约束 | 内容 |
|------|------|
| G1 覆盖 | $`\bigcup_i V_i=\Omega`$ |
| G2 缝相容 | 分片拼接时缝上归属一致 |
| G3 不过定 | 钉点总数与自由站点数相容 |
| $`\mathcal B`$ | 各指定子弧归属正确；渗漏等按故事验收 |

**输出（宽题三分类）**：无解 / 唯一解 / 多解（含正维解族）。宽题**不要求**先验给定「矩形/扇形/L」名称，也**不要求**全局初等闭式 $`p_i=p_i(\Omega)`$。

### 4.2　特设 $\Sigma$（收窄记录）

| 编号 | 特设 | 与宽题 W 的关系 |
|------|------|-----------------|
| $`\Sigma_2`$ | 路线 A：钉 $`M\in\Gamma_{12}`$，$`\lambda`$ 由 (1.3) 定 | 单钉外特征子题 |
| $`\Sigma_5`$ | 强规范 $`\mathcal B^{\ast}`$：渗漏边无 $`V_2`$ 严格内点 | 推论 M 的规范 |
| $`\Sigma_6`$ | $`n_{\min}`$ 相对 $`\mathcal B^{\ast}`$ 定义 | 定理 M★ |
| $`\Sigma_9`$ | 扇形 (I)(II)(III) 故事 | §8 |
| $`\Sigma_{10}`$ | 乘性泰森 $`\varphi_i=d/r_i`$ | 主模型 |

### 4.3　$`\mathfrak{Sol}^{\mathrm{mul}}`$

胞间分界线 $`\{\varphi_i=\varphi_j\}`$ 为阿波罗尼奥斯圆（$`r_i=r_j`$ 时退化为直线）。$`\partial\Omega`$ 的直边、光滑弧为域裁剪，一般不与某 $`\Gamma_{ij}`$ 重合。

$n=2$ 时边界审计用性价比差 $`\chi_{12}(x)=\varphi_2(x)-\varphi_1(x)`$。内禀违约集（无需人工标渗漏边）：

$$
\mathrm{Viol}(\mathcal S;\mathcal B):=\{s\in\partial\Omega:\ \sigma_{\mathcal S}(s)\neq\tau_{\mathcal B}(s)\}.
$$

**指标**（半代数维数）：

$$
\mathfrak I(f_D,\gamma,p_1,\tau;\mathfrak g):=\dim_{\mathbb R}\,\mathfrak{Sol}^{\mathrm{mul}}.
$$

### 定理 W.10（宽题代数归约）

宽题 W 的「存在配置 $`\mathcal S`$ 满足 $`\mathcal G\cup\mathcal B`$」等价于 $`\mathfrak{Sol}^W\neq\varnothing`$，其中 $`\mathfrak{Sol}^W`$ 为对有限 $N$ 与有限拓扑分支 $`\mathfrak g`$ 上 $`\mathfrak{Sol}^{\mathrm{mul,full}}`$ 的并集投影。

### 定理 W（宽题三分类）

| W 输出 | 解析判定 |
|--------|----------|
| 无解 | $`\mathfrak{Sol}^W=\varnothing`$ |
| 唯一 | $`\mathfrak{Sol}^W\neq\varnothing`$，本质参数仅一个等价类 |
| 多解 | 两解不等价，或某分支 $`\mathfrak I>0`$ |

### 定理 M★（同定理 1.1）

在路线 A、强规范 $`\mathcal B^{\ast}`$ 下：

$$
n_{\min}=\begin{cases}
2, & m_{\mathrm{eff}}=0,\\
3, & m_{\mathrm{eff}}\ge 1.
\end{cases}
$$

$`m_{\mathrm{eff}}\ge 1`$ 时第三面旗可取秩 3 配置 $`(p_3,r_3)=(p_2,\lambda)`$；引理 M★.0：此时**全局无** $`V_2`$ 严格内点，故所有渗漏边上的旗 2 开段同时消失。几何分离路线在 $`m_{\mathrm{eff}}\ge 2`$ 时可能需要 $`n_{\min}^{\mathrm{geo}}=2+m_{\mathrm{eff}}`$（推论 M.III–IV），与秩 3 路线并列。

### 推论 M 摘要（特设 $\Sigma$ 内可证片段）

| 段 | 内容 |
|----|------|
| M.I | $n=2$ 可行 $\Leftrightarrow$ $`\forall e\in\mathcal E,\ \min\chi_e\ge0`$；成立则 $`n_{\min}=2`$ |
| M.II | M.I 失败 $`\Rightarrow`$ 秩 3 共点 $`\Rightarrow`$ $`n_{\min}=3`$（步骤 3 恒成功） |
| M.III | 几何分离下界 $`n_{\min}^{\mathrm{geo}}\ge 2+m_{\mathrm{eff}}`$ |
| M.IV | 逐边守卫参数存在时 $`n_{\min}^{\mathrm{geo}}=2+m_{\mathrm{eff}}`$ |

### 4.10　已证与未证

| 已证 | 未证 |
|------|------|
| W-存在/唯一/多解归约（W.10 + W） | 弱 $`\mathcal B`$ 下 $`n_{\min}`$ |
| $`\mathrm{Viol}`$ 内禀化 | 一般 $`\Omega`$ 的 $`Z(\Omega)`$ 统一指标 |
| $`\Sigma_9`$ 扇形三分 | $`b_\eta`$ 全局初等闭式 |
| 定理 M★：$`n_{\min}\in\{2,3\}`$（秩 3） | 几何分离路线对一切 $`\Omega`$ 的存在性 |

### 4.4　假设 H

| 编号 | 内容 |
|------|------|
| H1 | $`\Omega`$ 有界，$`\partial\Omega=\bigcup_i\Gamma_i`$，每段局部 $`f_i(x,y)=0`$，$`\nabla f_i\neq0`$ |
| H2 | $`\mathcal E\subset\{\Gamma_i\}`$ 有限；每 $`e\in\mathcal E`$ 有 $C^1$ 参数化 $`\gamma_e`$ |
| H3 | 路线 A，$`F(x,y)=\frac{d(p_2)}{\lambda}-d(p_1)`$ |
| H4 | 单组 $`\mathcal S`$；拼接时 G2 |

### 4.5　方程系统 $`\mathfrak E^{\mathrm{mul}}`$

参数化 $`\partial D=\{\gamma(s):s\in\mathcal S\}`$。未知量 $`\mathcal U=(p_2,\ldots,p_N,r_1,\ldots,r_N)`$。

**(A) 外边界归属**　对 $`s\in\mathcal S`$，目标标签 $`\tau(s)`$：

$$
\varphi_{\tau(s)}(\gamma(s))\le \varphi_j(\gamma(s)),\quad \forall j\neq\tau(s).
$$

等号仅在 tie 集；平方去根号后为多项式等式/不等式。

**(B) 内胞邻接**　对每个预期相邻胞对 $(i,j)$，存在 $`x_{ij}\in D`$ 使 $`\varphi_i(x_{ij})=\varphi_j(x_{ij})`$，即

$$
\|x-p_i\|^2 r_j^2=\|x-p_j\|^2 r_i^2.
$$

**(C) 内顶点**　三条及以上 $`\Gamma_{ij}`$ 在 $D$ 内交于 $`x_\ast`$ 时，联立多个二次方程。

**(D) 可行性**　$`r_i>0`$；G3 钉点 $\Pi$；非退化 $`\mathcal U\notin\mathfrak D`$。

记 $`\mathcal U`$ 满足 $`\mathfrak E^{\mathrm{mul}}`$ 且 $`\mathrm{Viol}=\varnothing`$、$`\bigcup_i V_i=D`$，则

$$
\mathfrak{Sol}^{\mathrm{mul}}(f_D,\gamma,p_1,\tau;\mathfrak g):=\{\mathcal U\}.
$$

### 4.6　引理 M★.0

设 $`(p_3,r_3)=(p_2,\lambda)`$。则对一切 $`X\in\Omega`$，$X$ 不能为 $`V_2`$ 严格内点。

证明. $`X\in V_2`$ 严格内点当且仅当 $`d(X,p_2)/\lambda<d(X,p_1)`$ 且 $`d(X,p_2)/\lambda<d(X,p_3)/r_3`$。由 $`p_3=p_2`$、$`r_3=\lambda`$，后式化为 $`d(X,p_2)<d(X,p_2)`$，恒假。证毕.

### 4.7　定理 M★ 的证明

记 $`m_{\mathrm{eff}}=\lvert\{e\in\mathcal E:\ \mathcal S_2(e)\neq\varnothing\}\rvert`$。

若 $`m_{\mathrm{eff}}=0`$，则 $`\forall e\in\mathcal E`$，$`\min\chi_e\ge0`$，由 M.I 得 $`n_{\min}=2`$。

若 $`m_{\mathrm{eff}}\ge1`$，存在渗漏边，两面旗不满足 $`\mathcal B^{\ast}`$，故 $`n_{\min}\ge3`$。取 $`(p_3,r_3)=(p_2,\lambda)`$，由引理 M★.0，$`\mathcal E`$ 上无旗 2 严格开段；钉点 $M$ 仍在 $`\Gamma_{12}`$ 上。故 $n=3$ 可行，$`n_{\min}\le3`$。合并得 $`n_{\min}=3`$。证毕.

### 4.8　推论 M（完整）

**M.I**　强 $`\mathcal B^{\ast}`$ 下 $n=2$ 可行当且仅当 $`\forall e\in\mathcal E`$，$`\min_{t}\chi_e(t)\ge0`$，其中 $`\chi_e=F\circ\gamma_e`$。成立则 $`n_{\min}=2`$。

**M.II**　M.I 失败时，秩 3 共点 $`(p_3,r_3)=(p_2,\lambda)`$ 在路线 A 下实现强 $`\mathcal B^{\ast}`$，故 $`n_{\min}=3`$。

**M.III**　几何分离（$`p_3\not\equiv p_2`$）且守卫功能不相容时，$`n_{\min}^{\mathrm{geo}}\ge 2+m_{\mathrm{eff}}`$。

**M.IV**　若每条 $`e_i\in\mathcal E`$ 的 14.4 型守卫 $`(p_{2+i},r_{2+i})`$ 满足半径下界、$M$-安全且无交叉抢边，则 $`n_{\min}^{\mathrm{geo}}=2+m_{\mathrm{eff}}`$。

**M.V（障碍）**

| 障碍 | 条件 | 结论 |
|------|------|------|
| 非共点单站 | 对称 $`\mathcal B^{\ast}`$ + 单 $`p_3\not\equiv p_2`$ | 水平守卫不等式不可行 |
| 孔弧轴上分离 | $`p_3`$ 在轴上且 $`\chi(0)<0`$ | 几何分离无解 |
| 单钉吞不尽 | 单站钉 $`Q_e`$ | $`\mathcal S_2(e)`$ 内仍有 $`V_2`$ 严格点 |

### 4.9　$`Z(\Omega)`$（单连通）

设 $`\Omega`$ 单连通、边界分段 $C^1$。取 Riemann 映射 $`\Phi:\Omega\to\mathbb D`$。定义

$$
Z(\Omega):=\Bigl(\Phi,\ \theta\in\mathbb R/2\pi\mathbb Z\mapsto\gamma_{\mathrm{out}}(\theta):=\Phi^{-1}(e^{i\theta})\Bigr).
$$

$H\ge1$ 时 $`Z(\Omega):=(\Phi_{\mathrm{mc}},\{\theta_{\mathrm{out}}\},\{\theta_{\mathrm{in}}^{(h)}\})`$。扇形、矩形、L 形均为向 $`Z(\Omega)`$ 代入显式 $\Phi$ 的算例。

---

## 5　路线 A

### §3.1 参数化

取两旗、共线对称：

$$
p_1=(a,y_0),\quad r_1=1;\qquad p_2=(b,y_0),\quad r_2=\lambda>0. \qquad (1.2)
$$

未知数为 $`(b,\lambda)`$；**本征自由度** $`\mathrm{dof}=2`$。

### §3.2 钉点方程

边界点 $`M\in\partial\Omega`$ **钉**在两面旗分界线上，记 $`M\in\Gamma_{12}`$：

$$
\boxed{\frac{d(M,p_2)}{\lambda}=\frac{d(M,p_1)}{1}} \qquad (1.3)
$$

代数形式（与阿波罗尼奥斯圆等价）：

$$
\Psi_{12}(M):=\|M-p_1\|^2\lambda^2-\|M-p_2\|^2=0. \qquad (2.1)
$$

每个独立钉点增加 **1 个** 多项式等式 $E\gets E+1$。由维数引理（引理 2.1），若 $`E<\mathrm{dof}`$，解空间**必为正维**（无穷多解）。

### §3.3 路线 A 在宽题中的位置

路线 A **不是**宽题 W 的完整解，而是特设 $`\Sigma_2`$：先在外特征点 $M$ 钉 $`\Gamma_{12}`$，由 (1.3) 将 $`\lambda`$ 与 $b$ 联立。其充分性在扇形、矩形等算例已验证；与 v1「双钉 $A+M$」对照，见第八部分。

---

## 6　渗漏审计与 $`n_{\min}`$

### §4.1 性价比差 $`\chi`$

$N=2$ 时在边 $e$ 上比较

$$
\chi_{12}(x):=\varphi_2(x)-\varphi_1(x)=\frac{d(x,p_2)}{\lambda}-d(x,p_1).
$$

若存在开子段 $e'\subset e$ 使 $`\chi_{12}<0`$（旗 2 **严格**占优），而 $`\tau(e)=1`$（该边应归旗 1），则称 $e$ **渗漏**。

记渗漏边集 $`\mathcal E`$：

$$
m_{\mathrm{eff}}:=|\mathcal E|.
$$

$`m_{\mathrm{eff}}`$ 只依赖边界数据与候选 $`(p_2,\lambda)`$，**不依赖**「矩形/L/扇形」名称。

### §4.2 强规范 $`\mathcal B^{\ast}`$

**强规范** $`\mathcal B^{\ast}\subset\mathcal B`$：要求 $`\mathcal E=\varnothing`$，即渗漏边上无旗 2 严格开段。弱规范 $`\mathcal B`$ 只要求覆盖与 $`\tau`$ 一致，允许渗漏；推论 M 的 $`n_{\min}`$ 相对于 $`\mathcal B^{\ast}`$ 定义。

当 $`m_{\mathrm{eff}}\ge 1`$ 时，两旗无法在 $`\mathcal B^{\ast}`$ 下可行，必须增旗。

---

## 7　定点自由度

### §5.1 引理（维数下界，纯计数）

固定 $N=2$、路线 A 共线，$`\mathrm{dof}=2`$。设独立等式数（钉 + 额外代数约束）为 $E$。在一般位置钉方程独立时，

$$
\dim\mathfrak{Sol}\ge\max(0,\ \mathrm{dof}-E).
$$

| 情形 | 结论 |
|------|------|
| $`E<\mathrm{dof}`$ | **必**无穷解（至少 $`\mathrm{dof}-E`$ 维族） |
| $`E=\mathrm{dof}`$ | 可能有限解；需验相容性与 Viol |
| $`E>\mathrm{dof}`$ | 过定；一般无解或仅退化根 |

**引理 2.1 的证明**　$`\mathfrak{Sol}`$ 为半代数集；$E$ 个独立多项式方程至多降维 $E$。若 $`E<\mathrm{dof}`$，维数严格为正，不可能仅有限解。证毕.

### §5.2 定点个数速查表（$N=2$，共线，$`\mathrm{dof}=2`$）

| $`\lvert\Pi\rvert`$ | 额外约束 | $E$ vs dof | 解的形态 |
|----------|----------|--------------|----------|
| 0 | — | $0<2$ | **2 维族** |
| 1 | — | $1<2$ | **1 维族**（定理 1.2） |
| 1 | 面积份额 $`\eta`$ | $2=2$ | **0 维**；扇形上唯一 $`b_\eta`$（定理 1.3） |
| 1 | 固定 $`r_2`$ | $2=2$ | 有限（一般 0 或有限多点） |
| 2 | — | $2=2$ 但相容性差 | **一般无解**；扇形 v1：$b\in\{a,R^2/a\}$ |
| 2 | 无共线（$`\mathrm{dof}=3`$） | $2<3$ | 仍 **1 维族** |

固定 $N=2$ 共线：$`\lvert\Pi\rvert<\mathrm{dof}\Rightarrow\dim\mathfrak{Sol}\ge1`$；$`\lvert\Pi\rvert=\mathrm{dof}`$ 时查相容性；加独立第二方程（$`\eta`$ 等）可锁唯一。$`n_{\min}`$ 为最少旗数；$`\lvert\Pi\rvert`$ 为固定 $N$ 下的钉点个数。

---

## 8　定理 1.1

### 定理 1.1

在路线 A、强规范 $`\mathcal B^{\ast}`$ 下，

$$
n_{\min}=\begin{cases}
2, & m_{\mathrm{eff}}=0,\\
3, & m_{\mathrm{eff}}\ge 1.
\end{cases}
$$

当 $`m_{\mathrm{eff}}\ge 1`$ 时，第三面旗可取**秩 3 配置**

$$
(p_3,r_3)=(p_2,\lambda).
$$

### 证明

**下界**　$`m_{\mathrm{eff}}\ge 1`$ 时，存在边 $`e\in\mathcal E`$ 使 $`\tau(e)=1`$ 而 $`\chi_{12}<0`$ 于 $e$ 上开子段。两面旗无法消除该开段，故 $`n_{\min}\ge 3`$.

**上界（秩 3 构造）**　令 $`(p_3,r_3)=(p_2,\lambda)`$。$`X\in\Omega`$ 为旗 2 **严格**内点当且仅当

$$
\frac{d(X,p_2)}{\lambda}<\frac{d(X,p_1)}{1},
\qquad
\frac{d(X,p_2)}{\lambda}<\frac{d(X,p_3)}{\lambda}.
$$

由 $`p_3=p_2`$、$`r_3=\lambda`$，第二式化为 $`d(X,p_2)<d(X,p_2)`$，恒假。故**全局无旗 2 严格内点**，渗漏边上的旗 2 开段消失，$`\mathcal B^{\ast}`$ 成立。钉点 $M$ 仍满足 (1.3)。故 $n=3$ 可行，$`n_{\min}\le 3`$. 结合下界，$`m_{\mathrm{eff}}\ge 1\Rightarrow n_{\min}=3`$.

$`m_{\mathrm{eff}}=0`$　无渗漏，两旗已够，$`n_{\min}=2`$. 证毕.

### 注（秩 3 共点构造）

增旗的直接需求是在渗漏边 $e$ 上压制旗 2 的严格占优。几何分离路线（不同位置的 $`p_3`$）在 $`m_{\mathrm{eff}}\ge 2`$ 时可能需要 $`n_{\min}^{\mathrm{geo}}=4`$（顶/底守卫不相容）。秩 3 路线选择**同点复制** $`(p_2,\lambda)`$，使第三列比较退化，从谓词层面取消「旗 2 严格内点」，而钉点 $M$ 不移动。这是 $`\mathcal B^{\ast}`$ 下达到 $n=3$ 的**最小**代数增旗方式。扇形上 $`V_1\cup V_2=\Omega`$ 已由引理 2.3、5.1(iv) 成立；增旗机制来自矩形/L 边界的渗漏审计，而非补未覆盖区域。

### 算例审计（$`m_{\mathrm{eff}}\Rightarrow n_{\min}`$）

| 算例 | $`m_{\mathrm{eff}}`$ | $`n_{\min}`$（秩3） |
|------|----------------------|---------------------|
| 实心扇形 $OA,OB$ | 0 | 2 |
| 同心孔扇形 | 0 | 2 |
| 矩形仅顶边渗漏 | 1 | 3 |
| 矩形顶+底 | 2 | 3（秩3）/ 4（几何分离） |
| L 形边界 $`\mathcal B^{\ast}`$ | 2 | 3 / 4 |
| L 形竖直外臂验收 | 0 | 2 |

脚本 `nmin_decision_audit.py` 复验上表。

---

## 9　定理 1.2

### 定理 1.2

矩形 $`\Omega=[0,W]\times[0,H]`$，

$$
p_1=\Bigl(a,\frac{H}{2}\Bigr),\quad
p_2=\Bigl(b,\frac{H}{2}\Bigr),\quad
M=\Bigl(W,\frac{H}{2}\Bigr)\in\Gamma_{12}.
$$

则

$$
\boxed{\,b=W-(W-a)\,r_2\,},\qquad p_{2y}=\frac{H}{2},\quad r_2>0 \qquad (1.4)
$$

为解集的 **1 维族**；$`r_2`$（即 $`\lambda`$）为自由参数。

### 证明

在 $y=H/2$ 上 (1.3) 化为

$$
\frac{W-b}{r_2}=W-a
\quad\Longrightarrow\quad
b=W-(W-a)r_2.
$$

方程 1 个，未知 $`(b,r_2)`$ 两个，$`\mathrm{dof}-E=2-1=1`$，故 $`\dim\mathfrak{Sol}\ge 1`$。对任意 $`r_2>0`$ 取 $b$ 如上，(1.3) 成立。证毕.

### 数值例（$W=2,H=1,a=2/5$）

$$
p_{2x}=2-\frac{8}{5}r_2,\quad p_{2y}=\frac12,\quad r_2>0.
$$

`solve_land.py` 输出两条参数化分支（$`p_{2x}`$ 关于 $`r_2`$ 的闭式），与上式一致。再加面积份额 $`\mathrm{Area}(V_2)=\eta\cdot 2`$ 得第二个方程，$`\eta`$ 给定后 $b$ 唯一（扇形型单调性，定理 1.3 同型）。

### 双钉矩形

$W=2,H=1,a=2/5$ 时同时钉 $M=(2,1/2)$ 与 $(1,0)$，$`\lvert\Pi\rvert=2`$，$`E=2=\mathrm{dof}`$。在 $y=1/2$ 与 $y=0$ 上 (1.3) 分别为

$$
\frac{1-b}{r_2}=2-\frac25,\qquad
\frac{\frac12-0}{r_2}=\sqrt{\Bigl(\frac12-\frac25\Bigr)^2+\Bigl(\frac12\Bigr)^2}.
$$

左式给出 $`b=2-(8/5)r_2`$；代入右式得关于 $`r_2`$ 的相容方程。一般无非退化实解，与 $\Sigma_9$-I 同型。

### 定理 4.3（有限解的必要条件）

若 $`\dim\mathfrak{Sol}=0`$，则 $`E\ge\mathrm{dof}`$。若 $`E=\mathrm{dof}`$ 且钉方程组雅可比列满秩，则局部 $`\dim\mathfrak{Sol}=0`$。

证明. 必要式由引理 2.1 逆否；充分式由隐函数定理。证毕.

---

## 10　扇形与 $\Sigma_9$

### §8.1 扇形设定

$$
\Omega=\left\{(\rho,\varphi):\ 0\le \rho\le R,\ \lvert\varphi\rvert\le\alpha\right\},\qquad \alpha\in(0,\pi),\ R>0.
$$

边界：外弧 $`\Gamma_{\mathrm{out}}`$（$`\rho=R`$），径向边 $OA,OB$（$`\varphi=\pm\alpha`$）。角点

$$
A=(R\cos\alpha,\ R\sin\alpha).
$$

已知 $`p_1=(a,0)`$，$0<a<R$，$`r_1=1`$。对称 Ansatz $`p_2=(b,0)`$，$`r_2=\lambda`$。势力区

$$
V_2=\left\{x\in\Omega:\ \frac{d(x,p_2)}{\lambda}<d(x,p_1)\right\},\qquad V_1=\Omega\setminus V_2.
$$

路线 A：仅钉 $M=(R,0)$，$a<b<R$，

$$
\lambda=\frac{R-b}{R-a}.
$$

### §8.2 路线 A 结构（引理 5.1–5.2）

**引理 5.1**　(i) $`A\in V_1`$；(ii) 径向边无渗漏；(iii) $`V_1,V_2`$ 连通；(iv) $`V_1\cup V_2=\Omega`$。

证明. (i) 令 $P(\rho)=(\rho\cos\alpha,\rho\sin\alpha)$，

$$
f(\rho):=\|P(\rho)-p_1\|^2-\lambda^{-2}\|P(\rho)-p_2\|^2.
$$

在 $`\lambda=(R-b)/(R-a)`$ 下，$P(R)=A$。记 $`h=\|A-p_1\|=\sqrt{R^2+a^2-2aR\cos\alpha}`$。要证 $f(R)<0$，即 $`d(A,p_1)<d(A,p_2)/\lambda`$。代入 $`\lambda=(R-b)/(R-a)`$，等价于

$$
h(R-a)<(R-b)\sqrt{R^2+b^2-2bR\cos\alpha}.
$$

在 $a<b<R$、$\alpha\in(0,\pi)$ 下成立（路线 A 标准估计）。故 $`A\in V_1`$.

(ii) $f(\rho)$ 为 $\rho$ 的二次函数。$`f(0)=a^2-\lambda^{-2}b^2<0`$（$b>a$），且 $f$ 在 $(0,R)$ 无零点（否则 $OA$ 上出现 $`\Gamma_1`$ 内交，与引理 5.2 矛盾）。故径向边无渗漏。

(iii) $O=(0,0)$：$`d(O,p_1)=a`$，$`d(O,p_2)/\lambda=b(R-a)/(R-b)>a`$。$`V_1`$ 含 $O$ 与径向开射线，连通。

(iv) 引理 2.3。证毕.

**引理 5.2**　$`\Gamma_1`$ 与 $\rho=R$ 在 $M$ 内切。

证明. (1.3) 给出 $`M\in\Gamma_1`$；圆心在 $x$ 轴，与外圆在 $M$ 曲率匹配。证毕.

| 编号 | 性质 |
|------|------|
| S1 | $`M\in\Gamma_1`$ |
| S2 | $`A\in V_1`$；径向开射线 $\subset V_1$ |
| S3 | $`\Gamma_1`$ 与外弧在 $M$ **内切** |
| S4 | $`V_1,V_2`$ 连通 |
| S5 | $`V_1\cup V_2=\Omega`$ |
| S6 | $`A(b)=\mathrm{Area}(V_2)`$ 在 $(a,R)$ 连续、**严格单调减** |

脚本 `sector_sigma9_verify.py` 复验 S1–S6。

### §8.3 定理 1.3（面积份额 $`\Rightarrow`$ 唯一）

在路线 A 的 1 维族上，加

$$
\mathrm{Area}(V_2)=\eta\,\lvert\Omega\rvert,\qquad \eta\in(0,1),
$$

则存在唯一 $`b_\eta\in(a,R)`$。

**引理 5.3**　$`a<b_1<b_2<R\Rightarrow V_2(b_2)\cap\Omega\subseteq V_2(b_1)\cap\Omega`$。

证明. $`\psi_b(x)=\varphi_2(x)-\varphi_1(x)=0`$ 等价于

$$
\bigl((x_1-b)^2+x_2^2\bigr)(R-a)^2=\bigl((x_1-a)^2+x_2^2\bigr)(R-b)^2.
$$

展开得 $`(b-a)F_x(b)=0`$，$`F_x`$ 为一次式：

$$
\begin{aligned}
F_x(b)=&\ (R^2 a+R^2 b-2R^2 x_1-2Rab+2Rx_1^2+2Rx_2^2\\
&+2abx_1-ax_1^2-ax_2^2-bx_1^2-bx_2^2.
\end{aligned}
$$

故 $(a,R)$ 内至多一个 $b$ 使 $`\psi_b(x)=0`$。若 $`x\in V_2(b_2)`$ 且 $`\psi_{b_1}(x)\ge0`$，中间有零点，矛盾。证毕.

另一形式：

$$
b_\ast(x)=\frac{R^2 a-2R^2 x_1+2R x_1^2+2R x_2^2-a x_1^2-a x_2^2}
{R^2-2Ra+2ax_1-x_1^2-x_2^2}.
$$

**引理 5.4**　$A(b)$ 在 $(a,R)$ 连续、严格减；$`\lim_{b\uparrow R}A(b)=0`$；$`A_{\sup}:=\lim_{b\downarrow a}A(b)>0`$。

证明. 连续性由特征函数与控制收敛定理。$b\uparrow R$ 时 $`\lambda\to0^+`$，$`V_2\cap\Omega\to\{p_2\}`$，$A\to0$。$b\downarrow a$ 时月牙膨胀，$`A_{\sup}>0`$。单调性由引理 5.3。证毕.

**定理 1.3 的证明**　$A$ 值域 $`(0,A_{\sup})`$，连续严格减。对

$$
\eta\lvert\Omega\rvert\in(0,A_{\sup}),
$$

介值定理与严格单调给出唯一 $`b_\eta`$。证毕.

### §8.4 定理 $\Sigma_9$（扇形子故事穷举）

固定扇形、$`p_1`$、两旗乘性泰森。宽题输出**仅为**下表五行之一：

| 编号 | 子故事 | 硬约束 | W 输出 |
|------|--------|--------|--------|
| $\Sigma_9$-I | v1 双钉 | $`A,M\in\Gamma_1`$ | **无解** |
| $\Sigma_9$-II | 错故事 | $`M\in\Gamma_1`$ 且 $`A\in V_2`$ | **无解** |
| $\Sigma_9$-III | 路线 A | 仅 $`M\in\Gamma_1`$，无渗漏 | **多解**（1 维族） |
| $\Sigma_9$-IV | A + 份额 $`\eta`$ | 同 III + $`\mathrm{Area}(V_2)=\eta\lvert\Omega\rvert`$ | **唯一** |
| $\Sigma_9$-V | A + 第二钉 $A$ | 同 III + $`A\in\Gamma_1`$ | **无解** |

#### $\Sigma_9$-I 证明（v1 双钉无解）

(I) $`A\in\Gamma_1`$，(II) $`M\in\Gamma_1`$ 给出 $`\lambda=(R-b)/(R-a)`$。代入 (I) 并平方整理，

$$
p(b):=A_2b^2+A_1b+A_0=0,
$$

其中

$$
\begin{aligned}
A_2 &= 2aR(1-\cos\alpha),\\
A_1 &= -2R\Big[(R-a)\|A-p_1\| + (R-a)^2\cos\alpha\Big],\\
A_0 &= (R-a)^2 R^2 + \|A-p_1\|^2 R^2 - 2\|A-p_1\|^2 R(R-a).
\end{aligned}
$$

$`A_2>0`$（$\alpha\neq0$）。$b=a$ 时 $`\lambda=1`$，(I) 恒成立，为退化根。韦达定理给出另一根 $b=R^2/a$；由 $a<R$ 得 $R^2/a>R$。均不在 $(a,R)$。故 $`\mathfrak{Sol}=\varnothing`$.

#### $\Sigma_9$-II 证明

引理 5.1(i)：$a<b<R$ 时 $f(R)<0$，即 $`A\in V_1`$，与 $`A\in V_2`$ 矛盾。

#### $\Sigma_9$-III

任意 $b\in(a,R)$，$`\lambda(b)=(R-b)/(R-a)`$，引理 5.1–5.2 成立；$b$ 本质自由 $`\Rightarrow\dim=1`$.

#### $\Sigma_9$-IV

定理 1.3。

#### $\Sigma_9$-V

等价于在 III 上再加 (I)，同 $\Sigma_9$-I。

定理 $\Sigma_9$　扇形两旗乘性泰森下，输出为 $\Sigma_9$-I 至 V 之一；唯一仅 IV，多解仅 III。证毕.

### §8.5 v1 / v2 / 路线 A 对照

| 版本 | 钉点 | 角点 $A$ | $(a,R)$ 内解 |
|------|------|------------|----------------|
| v1 | $A+M$ | 在 $`\Gamma_1`$ | **无** |
| v2（错） | $M$ + 要 $`A\in V_2`$ | — | **无** |
| **路线 A** | 仅 $M$ | 在 $`V_1`$ | **1 维族**；$`\eta`$ 得唯一 |

---

## 11　矩形与 L 形

### §9.1 矩形渗漏

$`\Omega=[0,W]\times[0,H]`$，路线 A 单钉 $M=(W,H/2)$。

顶边 $`e_{\mathrm{top}}`$ 上 $`\chi_{12}`$ 在 $`Q=(x_Q,H)`$ 取极小；$W=2,H=1,a=0.5,b=0.7$ 时 $`x_Q\approx1.30`$，渗漏区间 $x\gtrsim(0.82,1.95)$。

旗 2 可严格占优而 $`\tau`$ 要求归旗 1 $`\Rightarrow m_{\mathrm{eff}}\ge 1`$（仅顶边审计时 $`m_{\mathrm{eff}}=1`$）。顶+底同时审计时 $`m_{\mathrm{eff}}=2`$，秩 3 路线 $`n_{\min}=3`$；几何分离路线 $`n_{\min}^{\mathrm{geo}}=4`$（顶/底守卫不相容）。

单钉闭式由定理 1.2、推论 4.1 给出。$`m_{\mathrm{eff}}=2`$ 时引入

$$
p_3=(\xi,H-\varepsilon),\quad p_4=(\xi,\varepsilon),\quad r_3=r_4=r
$$

压制顶底渗漏。矩形路线 C 证书：$\xi=9/10,\varepsilon=3/20,r=7/10$（`verify_guard_params.py`）。

### §9.2 L 形

$$
\Omega=([0,R_x]\times[0,w])\cup([0,w]\times[0,R_y]),\quad w\le R_x,R_y.
$$

标准参数 $`R_x=2,w=0.5,R_y=2`$，$`p_1=(0.4,0.25)`$，$`p_2=(0.6,0.25)`$，$`\lambda=7/8`$，只钉 $`M_x=(2,0.25)`$。

**定理 7.1**　只钉 $`M_x`$ 时：(1) $`V_1\cup V_2=\Omega`$；(2) $`M_x\in\Gamma_1`$；(3) $`M_y=(w/2,R_y)\in V_1`$（标准参数）；(4) 内凹角 $O=(0,0)$ 邻域 $\subset V_1$。

证明. (1) 引理 2.3。(2) (1.3)。(3) 数值：$`R_x=2,w=0.5,a_1=0.4,b_1=0.6`$ 时 $`d(M_y,p_1)<d(M_y,p_2)/\lambda`$ 于水平臂比较。(4) $`p_1`$ 最近。证毕.

| 现象 | 内容 |
|------|------|
| G1 | $`V_1\cup V_2=\Omega`$ |
| 内部 | 水平臂约 **73%** 面积为 $`V_2`$（面积故事，非边界 $`\chi`$） |
| 边界 $`\mathcal B^{\ast}`$ | 顶/底边 $x>w$ 段渗漏 $`\Rightarrow m_{\mathrm{eff}}=2`$，$`x_Q\approx0.92`$，渗漏区间 $\approx(0.55,1.95)$ |
| $`n_{\min}`$ | 秩 3：3；几何分离：4 |

**定理 7.4**　同时钉 $`M_x=(R_x,w/2)`$、$`M_y=(w/2,R_y)`$ 时，

$$
\lambda_x=\frac{R_x-b_1}{R_x-a_1},\qquad
\lambda_y=\frac{R_y-b_2}{R_y-a_2}.
$$

路线 A 仅一组 $`(p_2,\lambda)`$，故一般 $`\lambda_x\neq\lambda_y`$，$`\mathfrak{Sol}=\varnothing`$ 或退化为 $`p_2=p_1`$。

证明. 与 §8.4 $\Sigma_9$-I 双钉过定同构。证毕.

路线 C：$\xi=11/20,\varepsilon=1/10,r=9/10$（`verify_guard_params.py --shape L`）。

---

## 12　半代数

钉方程 $`\Psi_{ij}=0`$、$`r_i>0`$、边界条件 $`\chi_{12}\ge0`$ 于 $`\partial\Omega`$ 上均为多项式等式/不等式，故 $`\mathfrak{Sol}`$ 为半代数集；$`\dim\mathfrak{Sol}`$ 与分支数原则上可由 Tarski–Seidenberg 判定。

记 $`\mathfrak I=\dim_{\mathbb R}\mathfrak{Sol}^{\mathrm{mul}}`$，则

| $`\mathfrak I`$ | $`\lvert\mathfrak{Sol}\rvert`$ | 分类 |
|-------------|------------------|------|
| — | 0 | 无解 |
| 0 | 1 | 唯一 |
| 0 | $>1$ | 有限多解 |
| $\ge1$ | — | 无穷多解 |

`mul_tyson_solve.py` 对 $N=2$ 作 Gröbner 消元，输出有限解或参数化族（如定理 1.2），与引理 2.1 计数一致。

### 可判定性实验（凸域单钉）

| 实例 | $`\lvert\Pi\rvert`$ | 未知数 | 等式 | 剩余 dof | 提示 |
|------|----------|--------|------|----------|------|
| 矩形 + 单钉 | 1 | 3 | 1 | 2 | 多解 |
| 矩形 + 双钉 | 2 | 3 | 2 | 1 | 一般无解/退化 |
| 凸五边形 + 单钉 | 1 | 3 | 1 | 2 | 多解 |

（`decidability_convex_pin.py` 复现。）

---

## 13　算例

**例 1（矩形单钉）**　$W=2,H=1,a=2/5$：$`p_{2x}=2-(8/5)r_2`$，$`p_{2y}=1/2`$；$`r_2=0.08\Rightarrow p_{2x}=1.872`$；无空多解。

**例 2（扇形路线 A）**　$R=2,a=0.4,b=1.1$，$`\lambda=9/16`$；任意 $b\in(0.4,2)$ 合法。$`\eta=0.35\Rightarrow b_\eta\approx0.943`$。

**例 3（扇形 v1）**　根 $b=0.4,b=10$，无 $(0.4,2)$ 内解。

**例 4（L 形）**　$`p_2=(0.6,0.25)`$，$`\lambda=7/8`$；$`m_{\mathrm{eff}}=2\Rightarrow n_{\min}=3`$。

**例 5（矩形顶边审计）**　$W=2,H=1,a=0.5,b=0.7$，$`\lambda=(2-b)/(2-a)`$。顶边 $`\chi_{12}`$ 在 $`x_Q\approx1.30`$ 取极小；渗漏区间 $x\in(0.82,1.95)$；$`m_{\mathrm{eff}}=1`$。

**例 6（扇形面积锁）**　$R=2,a=0.4$，$`\eta=0.2`$ 时 $`b_\eta\approx0.72`$；$`\eta=0.5`$ 时 $`b_\eta\approx0.55`$（$A(b)$ 严格单调）。

**例 7（双钉矩形）**　$M=(2,1/2)$ 与 $(1,0)$ 联立 (1.3) 两式，无非退化解。

---

## 14　Laguerre 幂图与 MWVD

乘性泰森 $`\varphi_i=d/r_i`$ 的旗间分界线一般为阿波罗尼奥斯圆；Laguerre 幂图 $`\pi_i=\|x-p_i\|^2-r_i^2`$ 的旗间分界线为直线。本文主链采用乘性泰森。经典 MWVD 中加权 Voronoi 单元为连通开集；本文反问题在 $`\mathcal B^{\ast}`$ 下额外要求边界无渗漏，从而出现 $`m_{\mathrm{eff}}\Rightarrow n_{\min}`$ 与秩 3 构造，此非 MWVD 存在性理论本身的内容。

---

## 15　已证与未证

| 已证 | 未证 |
|------|------|
| $\Sigma_9$ 扇形三分 | 一般 $`\Omega`$ 的 $`Z(\Omega)`$ 统一指标 |
| $`m_{\mathrm{eff}}\Rightarrow n_{\min}\in\{2,3\}`$（秩 3） | 弱 $`\mathcal B`$ 下 $`n_{\min}`$ |
| dof vs $`\lvert\Pi\rvert`$ 解空间维数 | $`m_{\mathrm{eff}}\ge2`$ 几何分离完备性 |
| 矩形/L 标准算例与路线 C 证书 | $`b_\eta`$ 初等闭式 |

---

## 16　解空间分类

**定义 14.1**　输入 $`(\Omega,p_1,\tau)`$，输出 $`(p_2,\ldots,p_N,r_2,\ldots,r_N)`$ 及解空间类型：

| 类型 | 条件 |
|------|------|
| 无解 | $`\mathfrak{Sol}=\varnothing`$ 或 Viol 无可行实解 |
| 唯一 | $`\dim\mathfrak{Sol}=0`$ 且 $`\lvert\mathfrak{Sol}\cap\mathfrak D^c\rvert=1`$ |
| 有限多解 | $`\dim\mathfrak{Sol}=0`$，$`\lvert\mathfrak{Sol}\rvert>1`$，全离散 |
| 无穷多解 | $`\dim\mathfrak{Sol}\ge1`$，参数化闭式 |

**定义 14.2（解析解）**　边界匹配方程组的实根，或 Gröbner 消元给出的参数化闭式；允许正维解族。数值迭代不进证明链。

---

## 17　判定规格

| 步骤 | 条件 | 结论 |
|------|------|------|
| 1 | 给定 $`\Omega,p_1,\tau`$；路线 A | 未知 $`(b,\lambda)`$；钉点 $M$ 在 $`\Gamma_{12}`$ 上，得 (1.3) |
| 2 | 逐边 $`\chi_{12}`$ | $`m_{\mathrm{eff}}=\lvert\mathcal E\rvert`$ |
| 3 | $`m_{\mathrm{eff}}=0`$ | 定理 1.1：$`n_{\min}=2`$ |
| 4 | $`m_{\mathrm{eff}}\ge1`$ | 定理 1.1：秩 3，$`n_{\min}=3`$ |
| 5 | $`\lvert\Pi\rvert<\mathrm{dof}`$ | 引理 2.1：$`\dim\mathfrak{Sol}\ge1`$ |
| 6 | $`\lvert\Pi\rvert=\mathrm{dof}`$ | 查相容；双钉常无解（$\Sigma_9$-I） |
| 7 | 加 $`\eta`$ 或固定 $`r_2`$ | 定理 1.3 / $\Sigma_9$-IV：唯一 |

**情形对照**

| 情形 | 代数/审计事实 | 输出 |
|------|--------------|------|
| v1 双钉 $A+M$ | 引理 5.5：$`b\in\{a,R^2/a\}`$ | $\Sigma_9$-I，无解 |
| 单钉 $M$ | $`E=1<\mathrm{dof}=2`$ | 定理 1.2，1 维族 |
| 单钉 + $`\eta`$ | $`E=2=\mathrm{dof}`$，$A(b)$ 严格单调 | 定理 1.3，唯一 |
| 渗漏边存在 | $`m_{\mathrm{eff}}\ge1`$，两旗不满足 $`\mathcal B^{\ast}`$ | 定理 1.1，秩 3，$`n_{\min}=3`$ |
| 矩形/L/扇形 | 仅代入 $`\partial\Omega`$ 算 $`m_{\mathrm{eff}}`$ 与 dof | 不另立形状定理 |

---

## 18　路线 C

$`m_{\mathrm{eff}}\ge2`$ 时引入

$$
p_3=(\xi,H-\varepsilon),\quad p_4=(\xi,\varepsilon),\quad r_3=r_4=r.
$$

**M-安全**：$`d(M,p_1)\le d(M,p_k)/r_k`$。**渗漏边**：水平边上无 $`V_2`$ 严格内点。

| 形状 | $\xi$ | $\varepsilon$ | $r$ |
|------|---------|-----------------|-------|
| 矩形 | $9/10$ | $3/20$ | $7/10$ |
| L 形 | $11/20$ | $1/10$ | $9/10$ |

---

## 19　孔洞边界

设 $`\partial\Omega=\Gamma_{\mathrm{out}}\sqcup\bigsqcup_h\Gamma_{\mathrm{in}}^{(h)}`$。孔缘 $`\Gamma_{\mathrm{in}}^{(h)}`$ 的归属由 $`\tau`$ 指定（例如归 $`V_1`$）。带孔扇形上孔缘点 $N$ 与外弧点 $M$ 同时钉 $`\Gamma_1`$ 时，与 $\Sigma_9$-I 同型过定：联立两式 $`\lambda`$ 一般不相容，应只钉 $M$，孔缘归属由 (A) 推出。

---

## 20　算例推论索引

| 推论 | $`\Omega`$ | $`m_{\mathrm{eff}}`$ | $`n_{\min}`$ |
|------|----------|-------------------|------------|
| M.8 | 实心扇形 | 0 | 2 |
| M.9 | 同心孔扇形 | 0 | 2 |
| M.14 | 矩形顶底 | 2 | 3（秩3）/ 4（几何分离） |
| M.15 | L 形边界 | 2 | 3 / 4 |
| M.21 | 矩形/L 双水平边 | 2 | 3（秩3） |

---

## 21　算法 22.1（$\Sigma$ 内）

1. 构造 $`\mathcal E`$，计算各边 $`\chi_e`$ 的极小值。
2. 若 $`\forall e,\min\chi_e\ge0`$，输出 $`n_{\min}=2`$。
3. 否则取 $`(p_3,r_3)=(p_2,\lambda)`$，输出 $`n_{\min}=3`$。
4. 若拒绝秩 3，进入几何分离：下界 M.III，验证 M.IV 或终止于 M.V 障碍。
5. 固定 $N$，按 $`\lvert\Pi\rvert`$ 与 $`\mathrm{dof}`$ 分类 $`\mathfrak{Sol}`$；加 $`\eta`$ 锁唯一。
