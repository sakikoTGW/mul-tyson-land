# 乘性泰森地皮划分

> 给定平面闭域 \(\Omega\)、第一面旗 \(p_1\)、边界故事 \(\tau\)，反求其余旗位与半径，并判定解空间（无解 / 唯一 / 有限 / 无穷）。本文档给出主问题定义、推理链、核心定理与证明；仓库内脚本为数值复现，**非证明链**。

---

## 第一部分　几何对象与划分模型

### §1.1 地皮

地皮是一块平面**闭域** \(\Omega\subset\mathbb R^2\)：

- 外边界 \(\Gamma_{\mathrm{out}}\) 可由折线围成任意多边形（矩形、L 形、扇形、海星形等均为**代入同一框架的算例**，不另立形状定理）；
- 允许**孔洞**：\(\partial\Omega=\Gamma_{\mathrm{out}}\sqcup\bigsqcup_h\Gamma_{\mathrm{in}}^{(h)}\)，孔边界与外边界的归属均由 \(\tau\) 指定；
- 假设边界分段 \(C^1\)，除有限个角点外法向可定义。

### §1.2 旗与性价比

在平地上选定一点插**第一面旗** \(p_1\in\Omega\)，管辖半径取 \(r_1=1\)（归一化，不失一般性）。

尚需插入 \(N-1\) 面旗。第 \(i\) 面旗有位置 \(p_i\in\mathbb R^2\) 与管辖半径 \(r_i>0\)。配置记为

\[
\mathcal S=\{p_1,\ldots,p_N\}\cup\{r_1,\ldots,r_N\},\qquad r_1=1.
\]

对 \(x\in\mathbb R^2\)，到第 \(i\) 面旗的**性价比**（乘性泰森权）定义为

\[
\boxed{\ \varphi_i(x):=\frac{\|x-p_i\|}{r_i}\ }.
\]

\(x\) 归 \(\varphi_i\) **最小**的那面旗管辖。记势力区

\[
V_i:=\{x\in\Omega:\ \varphi_i(x)\le\varphi_j(x)\ \forall j\neq i\}
\]

（分界线上归属可任取一侧，零测度不影响面积与审计）。

**支配场**：

\[
\sigma_{\mathcal S}(x):=\operatorname*{argmin}_{1\le i\le N}\varphi_i(x).
\]

### §1.3 势力分界线：阿波罗尼奥斯圆

两旗 \(p_i,p_j\) 性价比相等时

\[
\frac{\|x-p_i\|}{r_i}=\frac{\|x-p_j\|}{r_j}
\quad\Longleftrightarrow\quad
\|x-p_i\|^2 r_j^2=\|x-p_j\|^2 r_i^2.
\]

这是平面上的**圆**（阿波罗尼奥斯圆）；仅当 \(r_i=r_j\) 时退化为垂直平分直线。因此乘性泰森胞腔边界**一般为圆弧**，与 Laguerre 幂图（等幂线为直线）不同——后者仅作对照，不进入主证明链。

### §1.4 边界故事 \(\tau\) 与可行性

**边界故事** \(\tau\)：对 \(\partial\Omega\) 的每条边（或子弧）\(e\)，指定 \(e\) 上除分界点外的目标归属 \(\tau(e)\in\{1,\ldots,N\}\)。

**程序性可行性**（整包约束 \(\mathcal G\) 之核心）：

| 编号 | 约束 | 形式 |
|------|------|------|
| G1 覆盖 | 整块地分完 | \(V_1\cup\cdots\cup V_N=\Omega\) |
| G2 缝相容 | 拼接片在缝上归属一致 | 分片 \(\Omega\) 时边界标签连续 |
| G3 不过定 | 钉点与自由未知数相容 | 见 §5 |
| \(\mathcal B\) | 边界故事成立 | \(\sigma_{\mathcal S}(s)=\tau(s)\) 对 \(s\in\partial\Omega\)（分界点除外） |

**违约集**（内禀审计，不必人工标渗漏边）：

\[
\mathrm{Viol}(\mathcal S;\tau):=\{s\in\partial\Omega:\ \sigma_{\mathcal S}(s)\neq\tau(s)\}.
\]

---

## 第二部分　反问题（主问题封口规格）

### §2.1 输入与输出

**输入**

1. 地皮 \(\Omega\)（顶点序列 + 可选孔洞）；
2. 第一面旗 \(p_1\in\Omega\)，\(r_1=1\)；
3. 边界故事 \(\tau\)（每条边归哪面旗）。

**过程**

4. 自动判定最少旗数 \(n_{\min}\)（定理 1，§4）；
5. 在固定 \(N\) 下建立方程组并消元。

**输出**

| 分类 | 数学含义 |
|------|----------|
| **无解** | \(\mathfrak{Sol}=\varnothing\) 或无可行实解满足 Viol |
| **唯一** | \(\dim\mathfrak{Sol}=0\) 且 \(|\mathfrak{Sol}\cap\mathfrak D^c|=1\) |
| **有限多解** | \(\dim\mathfrak{Sol}=0\) 且 \(|\mathfrak{Sol}|>1\)，全列 |
| **无穷多解** | \(\dim\mathfrak{Sol}\ge 1\)，给出参数化闭式与自由参 |

其中 \(\mathfrak D\) 为退化簇（\(p_i=p_j\)、\(r_i=0\) 等）。

### §2.2 解析解的含义

**解析解**指边界匹配代数方程组的**实根**，或 Gröbner 消元给出的**参数化闭式**；允许解族为正维，不要求全局初等闭式。数值迭代不进入证明链；仓库脚本用于复验。

---

## 第三部分　路线 A（\(N=2\)，共线 Ansatz）

### §3.1 参数化

取两旗、共线对称：

\[
p_1=(a,y_0),\quad r_1=1;\qquad p_2=(b,y_0),\quad r_2=\lambda>0.
\]

未知数为 \((b,\lambda)\)，**本征自由度** \(\mathrm{dof}=2\)。

### §3.2 钉点方程

边界点 \(M\in\partial\Omega\) **钉**在两面旗分界线上，记 \(M\in\Gamma_{12}\)：

\[
\boxed{\ \frac{d(M,p_2)}{\lambda}=\frac{d(M,p_1)}{1}\ }.
\tag{†}
\]

代数形式（与阿波罗尼奥斯圆等价）：

\[
\Psi_{12}(M):=\|M-p_1\|^2\lambda^2-\|M-p_2\|^2\cdot 1^2=0.
\]

每个独立钉点增加 **1 个** 多项式等式 \(E\gets E+1\)。由维数引理（§5.1），若 \(E<\mathrm{dof}\)，解空间**必为正维**（无穷多解）。

### §3.3 路线 A 的逻辑位置

路线 A **不是**宽题 W 的完整解，而是特设 \(\Sigma_2\)：先在外特征点 \(M\) 钉 \(\Gamma_{12}\)，由 (†) 将 \(\lambda\) 与 \(b\) 联立。其充分性在扇形、矩形等算例已验证；与 v1「双钉 \(A+M\)」对照，见 §7。

---

## 第四部分　渗漏审计与 \(n_{\min}\)

### §4.1 性价比差 \(\chi\)

\(N=2\) 时在边 \(e\) 上比较

\[
\chi_{12}(x):=\varphi_2(x)-\varphi_1(x)=\frac{d(x,p_2)}{\lambda}-d(x,p_1).
\]

若存在开子段 \(e'\subset e\) 使 \(\chi_{12}<0\)（旗 2 **严格**占优），而 \(\tau(e)=1\)（该边应归旗 1），则称 \(e\) **渗漏**。

记渗漏边集 \(\mathcal E\)，

\[
m_{\mathrm{eff}}:=|\mathcal E|.
\]

\(m_{\mathrm{eff}}\) 只依赖边界数据与候选 \((p_2,\lambda)\)，**不依赖**「矩形/L/扇形」名称。

### §4.2 强规范 \(\mathcal B^\*\)

**强规范** \(\mathcal B^\*\subset\mathcal B\)：要求 \(\mathcal E=\varnothing\)，即渗漏边上无旗 2 严格开段。弱规范 \(\mathcal B\) 只要求覆盖与 \(\tau\) 一致，允许渗漏；推论 M 的 \(n_{\min}\) 相对于 \(\mathcal B^\*\) 定义。

当 \(m_{\mathrm{eff}}\ge 1\) 时，两旗无法在 \(\mathcal B^\*\) 下可行，必须增旗。

---

## 第五部分　定点自由度与解空间维数

### §5.1 引理（维数下界，纯计数）

固定 \(N=2\)、路线 A 共线，\(\mathrm{dof}=2\)。设独立等式数（钉 + 额外代数约束）为 \(E\)。在一般位置钉方程独立时，

\[
\dim\mathfrak{Sol}\ge\max(0,\ \mathrm{dof}-E).
\]

| 情形 | 结论 |
|------|------|
| \(E<\mathrm{dof}\) | **必**无穷解（至少 \(\mathrm{dof}-E\) 维族） |
| \(E=\mathrm{dof}\) | 可能有限解；需验相容性与 Viol |
| \(E>\mathrm{dof}\) | 过定；一般无解或仅退化根 |

*证明*：\(\mathfrak{Sol}\) 为半代数集；\(E\) 个独立多项式方程至多降维 \(E\)。∎

### §5.2 定点个数速查表（\(N=2\)，共线，\(\mathrm{dof}=2\)）

| \(|\Pi|\) | 额外约束 | \(E\) vs dof | 解的形态 |
|----------|----------|--------------|----------|
| 0 | — | \(0<2\) | **2 维族** |
| 1 | — | \(1<2\) | **1 维族**（定理 2） |
| 1 | 面积份额 \(\eta\) | \(2=2\) | **0 维**；扇形上唯一 \(b_\eta\)（定理 3） |
| 1 | 固定 \(r_2\) | \(2=2\) | 有限（一般 0 或有限多点） |
| 2 | — | \(2=2\) 但相容性差 | **一般无解**；扇形 v1：\(b\in\{a,R^2/a\}\) |
| 2 | 无共线（dof\(=3\)） | \(2<3\) | 仍 **1 维族** |

**读图规则**

1. 数 \(\mathrm{dof}\)（是否共线 Ansatz）；
2. 数钉点 \(|\Pi|\) 与计划加的独立等式；
3. \(|\Pi|<\mathrm{dof}\Rightarrow\) 必无穷；
4. \(|\Pi|=\mathrm{dof}\Rightarrow\) 查相容性（双钉常失败）；
5. 若仍无穷，加第 \(\mathrm{dof}\) 个独立等式（\(\eta\)、固定半径等）\(\Rightarrow\) 唯一。

**与 \(n_{\min}\) 的分工**：\(n_{\min}\) 回答「最少几面旗」；本节回答「固定 \(N\) 下要钉几个点才有限/唯一」。二者正交。

---

## 第六部分　定理 1（最少旗数，秩 3）

### 定理 1

在路线 A、强规范 \(\mathcal B^\*\) 下，

\[
n_{\min}=\begin{cases}
2, & m_{\mathrm{eff}}=0,\\[6pt]
3, & m_{\mathrm{eff}}\ge 1.
\end{cases}
\]

当 \(m_{\mathrm{eff}}\ge 1\) 时，第三面旗取**秩 3 配置**

\[
(p_3,r_3)=(p_2,\lambda).
\]

### 证明

**下界**：\(m_{\mathrm{eff}}\ge 1\) 时，存在边 \(e\in\mathcal E\) 使 \(\tau(e)=1\) 而 \(\chi_{12}<0\) 于 \(e\) 上开子段。两面旗无法消除该开段，故 \(n_{\min}\ge 3\)。

**上界（秩 3 构造）**：令 \((p_3,r_3)=(p_2,\lambda)\)。\(X\in\Omega\) 为旗 2 **严格**内点当且仅当

\[
\frac{d(X,p_2)}{\lambda}<\frac{d(X,p_1)}{1},
\qquad
\frac{d(X,p_2)}{\lambda}<\frac{d(X,p_3)}{\lambda}.
\]

由 \(p_3=p_2\)、\(r_3=\lambda\)，第二式化为 \(d(X,p_2)<d(X,p_2)\)，恒假。故**全局无旗 2 严格内点**，渗漏边上的旗 2 开段消失，\(\mathcal B^\*\) 成立。钉点 \(M\in\Gamma_{12}\) 仍满足 (†)。故 \(n=3\) 可行，\(n_{\min}=3\)。

**\(m_{\mathrm{eff}}=0\)**：无渗漏，两旗已够，\(n_{\min}=2\)。∎

### 动机：为何共点第三旗

增旗的直接需求是在渗漏边 \(e\) 上压制旗 2 的严格占优。几何分离路线（不同位置的 \(p_3\)）在 \(m_{\mathrm{eff}}\ge 2\) 时可能需 \(n_{\min}^{\mathrm{geo}}=4\)（守卫不相容，§9）。秩 3 路线选择**同点复制** \((p_2,\lambda)\)，使第三列比较退化，从谓词层面取消「旗 2 严格内点」，而不移动 \(M\in\Gamma_{12}\)。这是 \(\mathcal B^\*\) 下达到 \(n=3\) 的**最小**代数增旗方式。

### 算例审计（\(m_{\mathrm{eff}}\) \(\Rightarrow\) \(n_{\min}\)）

| 算例 | \(m_{\mathrm{eff}}\) | \(n_{\min}\)（秩3） |
|------|----------------------|---------------------|
| 实心扇形 \(OA,OB\) | 0 | 2 |
| 同心孔扇形 | 0 | 2 |
| 矩形仅顶边渗漏 | 1 | 3 |
| 矩形顶+底 | 2 | 3（秩3）/ 4（几何分离） |
| L 形边界 \(\mathcal B^\*\) | 2 | 3 / 4 |
| L 形竖直外臂验收 | 0 | 2 |

脚本 `nmin_decision_audit.py` 复验上表。

---

## 第七部分　定理 2（矩形单钉：1 维族）

### 定理 2

矩形 \(\Omega=[0,W]\times[0,H]\)，

\[
p_1=\Bigl(a,\frac H2\Bigr),\quad
p_2=\Bigl(b,\frac H2\Bigr),\quad
M=\Bigl(W,\frac H2\Bigr)\in\Gamma_{12}.
\]

则

\[
\boxed{\,b=W-(W-a)\,r_2\,},\qquad p_{2y}=\frac H2,\quad r_2>0
\]

为解集的 **1 维族**；\(r_2\)（即 \(\lambda\)）为自由参数。

### 证明

在 \(y=H/2\) 上 (†) 化为

\[
\frac{W-b}{r_2}=W-a
\quad\Longrightarrow\quad
b=W-(W-a)r_2.
\]

方程 1 个，未知 \((b,r_2)\) 两个，\(\mathrm{dof}-E=2-1=1\)，故 \(\dim\mathfrak{Sol}\ge 1\)。对任意 \(r_2>0\) 取 \(b\) 如上， (†) 成立。∎

### 数值例（\(W=2,H=1,a=2/5\)）

\[
p_{2x}=2-\frac{8}{5}r_2,\quad p_{2y}=\frac12,\quad r_2>0.
\]

`solve_land.py` 输出两条参数化分支（\(p_{2x}\) 关于 \(r_2\) 的闭式），与上式一致。再加面积份额 \(\mathrm{Area}(V_2)=\eta\cdot 2\) 得第二个方程，\(\eta\) 给定后 \(b\) 唯一（扇形型单调性，定理 3 同型）。

### 双钉矩形

若同时钉 \(M=(W,H/2)\) 与底边中点 \((W/2,0)\)（\(|\Pi|=2,E=2=\mathrm{dof}\)），联立两 \(\Psi_{12}=0\) 与定理 9.5 / 扇形 v1 同型，**一般无非退化实解**。

---

## 第八部分　扇形与定理 Σ9（完整三分类）

### §8.1 扇形设定

\[
\Omega=\{(ρ,\varphi):\ 0\le ρ\le R,\ |\varphi|\le\alpha\},\quad \alpha\in(0,\pi),\ R>0.
\]

边界：外弧 \(\Gamma_{\mathrm{out}}\)（\(ρ=R\)），径向边 \(OA,OB\)（\(\varphi=\pm\alpha\)）。角点 \(A=(R\cos\alpha,R\sin\alpha)\)。

已知 \(p_1=(a,0)\)，\(0<a<R\)，\(r_1=1\)。对称 Ansatz \(p_2=(b,0)\)，\(r_2=\lambda\)。

\[
V_2=\{x\in\Omega:\ d(x,p_2)/\lambda<d(x,p_1)\},\quad V_1=\Omega\setminus V_2.
\]

### §8.2 路线 A 结构定理（摘要）

在仅钉 \(M=(R,0)\)、\(a<b<R\)、\(\lambda=(R-b)/(R-a)\) 时：

| 编号 | 性质 |
|------|------|
| S1 | \(M\in\Gamma_1\) |
| S2 | \(A\in V_1\)；径向开射线 \(\subset V_1\) |
| S3 | \(\Gamma_1\) 与外弧在 \(M\) **内切** |
| S4 | \(V_1,V_2\) 连通 |
| S5 | \(V_1\cup V_2=\Omega\)（两面旗已盖满，无空白） |
| S6 | \(A(b)=\mathrm{Area}(V_2)\) 在 \((a,R)\) 上连续、**严格单调减** |

*证明要点*：径向边函数 \(f(ρ)\) 在 \(ρ=R\) 处 \(f(R)<0\)（\(A\in V_1\)）；阿波罗尼奥斯圆与外圆内切；\(V_1\) 含原点路径连通。*脚本* `sector_sigma9_verify.py`。

### §8.3 定理 3（面积份额 \(\Rightarrow\) 唯一）

在路线 A 的 1 维族上，加

\[
\mathrm{Area}(V_2)=\eta\,|\Omega|,\qquad \eta\in(0,1),
\]

则存在唯一 \(b_\eta\in(a,R)\)。

### 证明

\(A(b):=\mathrm{Area}(V_2(b))\) 在 \((a,R)\) 上连续；\(b\uparrow R\) 时 \(\lambda\to0^+\)，\(A(b)\to0\)；\(b\downarrow a\) 时 \(A(b)\to A_{\sup}>0\)。引理：\(a<b_1<b_2<R\Rightarrow V_2(b_2)\cap\Omega\subseteq V_2(b_1)\cap\Omega\)（单零点引理，§8.10 代数版）。故 \(A\) **严格单调减**，值域 \((0,A_{\sup})\)。介值定理 \(\Rightarrow\ \exists b_\eta\)；严格单调 \(\Rightarrow\) 唯一。∎

### §8.4 定理 Σ9（扇形子故事穷举）

固定扇形、\(p_1\)、两旗乘性泰森。宽题输出**仅**为下表五行之一：

| 编号 | 子故事 | 硬约束 | W 输出 |
|------|--------|--------|--------|
| **Σ9-I** | v1 双钉 | \(A,M\in\Gamma_1\) | **无解** |
| **Σ9-II** | 错故事 | \(M\in\Gamma_1\) 且 \(A\in V_2\) | **无解** |
| **Σ9-III** | 路线 A | 仅 \(M\in\Gamma_1\)，无渗漏 | **多解**（1 维族） |
| **Σ9-IV** | A + 份额 \(\eta\) | 上 + \(\mathrm{Area}(V_2)=\eta|\Omega|\) | **唯一** |
| **Σ9-V** | A + 第二钉 \(A\) | 上 + \(A\in\Gamma_1\) | **无解** |

#### Σ9-I 证明（v1 双钉无解）

(I)(II) 联立化为 \(p(b)=A_2b^2+A_1b+A_0=0\)。韦达定理给出**恰两根** \(b\in\{a,\ R^2/a\}\)：\(b=a\) 退化（\(p_2=p_1\)）；\(b=R^2/a>R\)。均不在 \((a,R)\)。故 \(\mathfrak{Sol}=\varnothing\)。

#### Σ9-II 证明

定理 8.1(i)：\(a<b<R\) 时 \(f(R)<0\)，即 \(A\in V_1\)，与 \(A\in V_2\) 矛盾。

#### Σ9-III

任意 \(b\in(a,R)\)，\(\lambda(b)=(R-b)/(R-a)\)，定理 8.1–8.6 成立；\(b\) 本质自由 \(\Rightarrow\dim=1\)。

#### Σ9-IV

定理 3。

#### Σ9-V

等价于在 III 上再加 (I)，同 Σ9-I。

### §8.5 v1 / v2 / 路线 A 对照

| 版本 | 钉点 | 角点 \(A\) | \((a,R)\) 内解 |
|------|------|------------|----------------|
| v1 | \(A+M\) | 在 \(\Gamma_1\) | **无** |
| v2（错） | \(M\) + 要 \(A\in V_2\) | — | **无** |
| **路线 A** | 仅 \(M\) | 在 \(V_1\) | **1 维族**；\(\eta\) 得唯一 |

**推论**：扇形上两面旗已**全覆盖** \(\Omega\)（定理 8.6）；增第三旗的动机不是「补空白」，而是非扇形算例的 \(\mathcal B^\*\)（\(m_{\mathrm{eff}}\ge 1\)）或孔洞/非凸拓扑。

---

## 第九部分　矩形与 L 形（算例代入）

### §9.1 矩形渗漏

\(\Omega=[0,W]\times[0,H]\)，路线 A 单钉 \(M=(W,H/2)\)。顶边 \(e_{\mathrm{top}}\) 上旗 2 可严格占优而 \(\tau\) 要求归旗 1 \(\Rightarrow m_{\mathrm{eff}}\ge 1\)（仅顶边审计时 \(m_{\mathrm{eff}}=1\)）。顶+底同时审计时 \(m_{\mathrm{eff}}=2\)，秩 3 路线 \(n_{\min}=3\)；几何分离路线 \(n_{\min}^{\mathrm{geo}}=4\)（顶/底守卫不相容）。

路线 C 守卫参数（§21，非证明链）：矩形 \(\xi=9/10,\varepsilon=3/20,r=7/10\) 可复验 `verify_guard_params.py`。

### §9.2 L 形

\[
\Omega=([0,R_x]\times[0,w])\cup([0,w]\times[0,R_y]),\quad w\le R_x,R_y.
\]

标准参数 \(R_x=2,w=0.5,R_y=2\)，\(p_1=(0.4,0.25)\)，\(p_2=(0.6,0.25)\)，\(\lambda=7/8\)，只钉 \(M_x=(2,0.25)\)。

| 现象 | 内容 |
|------|------|
| G1 | \(V_1\cup V_2=\Omega\) |
| 内部 | 水平臂约 **73%** 面积为 \(V_2\)（面积故事，非边界 \(\chi\)） |
| 边界 \(\mathcal B^\*\) | 顶/底边 \(x>w\) 段渗漏 \(\Rightarrow m_{\mathrm{eff}}=2\) |
| \(n_{\min}\) | 秩3：3；几何分离：4 |

路线 C：\(\xi=11/20,\varepsilon=1/10,r=9/10\)（`verify_guard_params.py --shape L`）。

---

## 第十部分　推理主线（如何推到封口）

以下链条将「地皮插旗」反问题收束为可判定规格，各步均为前述定理的拼接，而非启发式补丁。

```
Ω, p₁, τ
    ↓  乘性泰森 φᵢ=d/rᵢ；Viol=σ≠τ
路线 A：p₂=(b,y₀), r₂=λ；钉 M∈Γ₁₂ → (†)
    ↓  逐边 χ₁₂；得 m_eff
m_eff=0 ──→ n_min=2（定理 1）
m_eff≥1 ──→ 秩3 (p₃,r₃)=(p₂,λ) → n_min=3（定理 1）
    ↓  固定 N；数 dof 与 |Π|
|Π|<dof ──→ 无穷解（引理 §5.1）
|Π|=dof ──→ 查相容（双钉常无解，Σ9-I）
+η 或 固定 r₂ ──→ 唯一（定理 3，Σ9-IV）
    ↓
输出：无解 / 唯一 / 有限 / 无穷 + 闭式
```

**关键转折**

1. **v1 \(\to\) 路线 A**：双钉 \(A+M\) 代数上仅给出 \(b\in\{a,R^2/a\}\)（Σ9-I），故改钉 \(M\)  alone，释放 1 维族。
2. **无穷 \(\to\) 唯一**：单钉仍 1 维（定理 2）；加**独立**第二方程（面积份额 \(\eta\)）使 \(E=\mathrm{dof}\)，扇形上单调性锁唯一（定理 3）。
3. **两旗 \(\to\) 三旗**：渗漏边存在时两旗不可 \(\mathcal B^\*\)；秩 3 共点增旗为最小修复（定理 1），与「盖满」无关（扇形已盖满）。
4. **形状名 \(\to\) 审计**：不分支「矩形定理/L 定理」，只代入 \(\Omega\) 边界算 \(m_{\mathrm{eff}}\) 与 \(\mathrm{dof}\)。

---

## 第十一部分　算法规格（工具层）

**输入**：顶点序列、\(p_1\)、\(\tau\)。

**步骤**

1. 建 `ProblemSpec`：边界 + 自动或手动钉点 + 内点见证；
2. Gröbner 消元（`mul_tyson_solve.py`）\(\Rightarrow\) 有限解 + 参数化族；
3. `audit_m_eff` \(\Rightarrow m_{\mathrm{eff}}\Rightarrow n_{\min}\)；
4. 若 \(m_{\mathrm{eff}}\ge 1\)，秩 3 验证 \((p_3,r_3)=(p_2,\lambda)\)；
5. 分类 \(\mathfrak{Sol}\)：Viol PASS 采样 + 维数计数；
6. 输出报告（`solve_land.py`）。

**退化处理**：参数族采样时优先共线、非贴边代表点（避免角点退化解误导作图）。

---

## 附录 A　脚本索引

| 脚本 | 对应数学环节 |
|------|----------------|
| `solve_land.py` | 封口：\(\Omega,p_1,\tau\to n_{\min}\) + 解空间 |
| `mul_tyson_solve.py` | \(\Psi_{ij}\) 消元 |
| `mul_tyson_viz.py` | 势力区可视化 |
| `sector_sigma9_verify.py` | Σ9-I–V 符号/数值 |
| `nmin_decision_audit.py` | 定理 1 算例表 |
| `verify_guard_params.py` | \(m_{\mathrm{eff}}\ge 2\) 路线 C |
| `decidability_convex_pin.py` | dof 实验 |
| `compare_boundary_story.py` | 与 CGAL/容量泰森对照 |
| `plot_sector_figures.py` | 扇形标准算例图 |
| `plot_l_shape_figures.py` | L 形标准算例图 |

```bash
pip install -r requirements.txt
python solve_land.py --vertices "0,0 2,0 2,1 0,1" --p1 "0.4,0.5" --tau "0:1 1:2 2:1 3:1"
python sector_sigma9_verify.py
python nmin_decision_audit.py
python verify_guard_params.py
```

---

## 附录 B　未纳入本仓库的内容

- **主文稿**（完整证明、共形工具 \(Z(\Omega)\)、阶段 III 一般 \(\Omega\)）：本地重写后单独发布；
- **示意图**（`figures/`、`*.png`）：本地生成，默认不推送。

---

## 参考文献（模型对照）

乘性加权 Voronoi / 阿波罗尼奥斯圆为经典 MWVD 结构；本文贡献在于 **\((\Omega,p_1,\tau)\) 反问题** 的解空间分类、\(m_{\mathrm{eff}}\Rightarrow n_{\min}\) 与定点 dof 计数的**封口规格**，而非重新发现圆的分界线形状。
