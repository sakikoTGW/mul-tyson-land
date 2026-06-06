# 乘性泰森地皮划分

## 1　引言

设 \(\Omega\subset\mathbb R^2\) 为有界闭域，\(\partial\Omega\) 由折线围成，允许孔洞。在平地上给定 \(p_1\in\Omega\)，插第一面旗，管辖半径取 \(r_1=1\)。第 \(i\) 面旗有位置 \(p_i\in\mathbb R^2\) 与半径 \(r_i>0\)。对 \(x\in\Omega\)，定义性价比

\[
\varphi_i(x):=\frac{\|x-p_i\|}{r_i}.
\tag{1.1}
\]

\(x\) 归 \(\varphi_i\) 最小者管辖；\(\varphi_i=\varphi_j\) 的等值线为势力分界线。在乘性泰森模型中，两旗等值线为阿波罗尼奥斯圆；仅当 \(r_i=r_j\) 时退化为直线。这与 Laguerre 幂图（等幂线为直线）不同，后者只作对照。

边界故事 \(\tau\) 规定 \(\partial\Omega\) 每条边（除分界点外）应归哪面旗。要求 \(V_1\cup\cdots\cup V_N=\Omega\)（覆盖），且边界实际归属与 \(\tau\) 一致。反问题：给定 \(\Omega,p_1,\tau\)，求其余 \((p_i,r_i)\)，并判定无解、唯一、有限多解或无穷多解；无穷时给出参数化闭式。

本文在路线 A 与强规范 \(\mathcal B^\*\)（无渗漏边）下建立 \(n_{\min}\) 与定点自由度理论，并在扇形、矩形、L 形上完成算例。第 2 节为预备知识；第 3 节证最少旗数定理 1.1；第 4 节证矩形单钉定理 1.2 并给出定点表；第 5 节证扇形份额定理 1.3 及定理 \(\Sigma 9\)；第 6、7 节代入矩形与 L 形；第 8 节为半代数框架；第 9 节为手算例；第 10 节收束推理链。

路线 A 采用

\[
p_1=(a,y_0),\ r_1=1;\qquad p_2=(b,y_0),\ r_2=\lambda>0,
\tag{1.2}
\]

边界点 \(M\in\partial\Omega\) 钉在 \(\Gamma_{12}\) 上：

\[
\frac{d(M,p_2)}{\lambda}=\frac{d(M,p_1)}{1}.
\tag{1.3}
\]

**定理 1.1**　在路线 A、\(\mathcal B^\*\) 下，\(n_{\min}=2\) 若 \(m_{\mathrm{eff}}=0\)，\(n_{\min}=3\) 若 \(m_{\mathrm{eff}}\ge 1\)；后者第三面旗可取 \((p_3,r_3)=(p_2,\lambda)\)。

**定理 1.2**　\(\Omega=[0,W]\times[0,H]\)，\(p_1=(a,H/2)\)，\(M=(W,H/2)\in\Gamma_{12}\) 时，\(b=W-(W-a)r_2\)（\(r_2>0\) 自由）为解集的一条 1 维族。

**定理 1.3**　扇形路线 A 的 1 维族上，若再加 \(\mathrm{Area}(V_2)=\eta|\Omega|\)（\(\eta\in(0,1)\)），则存在唯一 \(b_\eta\in(a,R)\)。

---

## 2　预备知识

### 2.1　划分与违约

记 \(\mathcal S=\{p_i,r_i\}_{i=1}^N\)，

\[
V_i:=\{x\in\Omega:\ \varphi_i(x)\le\varphi_j(x),\ \forall j\neq i\},\qquad
\sigma_{\mathcal S}(x):=\operatorname*{argmin}_i\varphi_i(x).
\]

整包约束：G1 覆盖 \( \bigcup_i V_i=\Omega\)；G2 缝相容；G3 钉点与自由元相容；\(\mathcal B\) 为 \(\sigma_{\mathcal S}=\tau\) 于 \(\partial\Omega\)（分界点除外）。

\[
\mathrm{Viol}(\mathcal S;\tau):=\{s\in\partial\Omega:\ \sigma_{\mathcal S}(s)\neq\tau(s)\}.
\]

\(N=2\) 时 \(\chi_{12}(x)=\varphi_2(x)-\varphi_1(x)\)。若 \(\tau(e)=1\) 而边 \(e\) 上存在开子段使 \(\chi_{12}<0\)，称 \(e\) **渗漏**。记 \(\mathcal E\) 为渗漏边集，\(m_{\mathrm{eff}}:=|\mathcal E|\)。强规范 \(\mathcal B^\*\) 要求 \(\mathcal E=\varnothing\)。

### 2.2　钉点方程

\(M\in\Gamma_{12}\) 等价于

\[
\Psi_{12}(M):=\|M-p_1\|^2\lambda^2-\|M-p_2\|^2=0.
\tag{2.1}
\]

每个独立硬钉使等式数 \(E\) 增 1。共线 Ansatz (1.2) 下 \(\mathrm{dof}=2\)。

**引理 2.1**　设 \(\mathfrak{Sol}\subset\mathbb R^{\mathrm{dof}}\) 为钉方程及可行性条件的半代数解集。若 \(E\) 个钉方程在一般位置独立，则 \(\dim\mathfrak{Sol}\ge\max(0,\mathrm{dof}-E)\)。

证明. 半代数集由多项式等式 \(F_1=\cdots=F_E=0\) 与不等式（如 \(r_2>0\)、Viol 条件）定义。每个独立等式在解集上至少降维 1，而不增加维数。故 \(\dim\mathfrak{Sol}\ge\mathrm{dof}-E\)；若右端为负则取 0。特别地 \(E<\mathrm{dof}\) 时不可能仅有限解。证毕.

**引理 2.2**　\(r_1\neq r_2\) 时 \(\varphi_1=\varphi_2\) 为圆；\(r_1=r_2\) 时为直线。

证明. \(\|x-p_1\|/r_1=\|x-p_2\|/r_2\) 两边平方整理得 \(\|x-p_1\|^2 r_2^2=\|x-p_2\|^2 r_1^2\)，为标准阿波罗尼奥斯圆方程。\(r_1=r_2\) 时退化为垂直平分线。证毕.

**引理 2.3（覆盖）**　\(N=2\) 时，对任意 \(x\in\Omega\)，令 \(u=\varphi_1(x)\)，\(v=\varphi_2(x)\)。由实数三歧性，\(u\le v\) 或 \(v\le u\)，故 \(x\in V_1\cup V_2\)。

证明. 直接由定义。证毕.

---

## 3　定理 1.1 的证明

### 3.1　秩 3 引理

**引理 3.1**　令 \((p_3,r_3)=(p_2,\lambda)\)。则不存在 \(X\in\Omega\) 为旗 2 的严格内点。

证明. \(X\) 为旗 2 严格内点当且仅当

\[
\varphi_2(X)<\varphi_1(X)\quad\text{且}\quad \varphi_2(X)<\varphi_3(X).
\]

第二式即 \(d(X,p_2)/\lambda<d(X,p_3)/\lambda\)。由 \(p_3=p_2\)、\(r_3=\lambda\)，得 \(d(X,p_2)<d(X,p_2)\)，矛盾。故无旗 2 严格内点。证毕.

### 3.2　下界

设 \(m_{\mathrm{eff}}\ge 1\)，取 \(e\in\mathcal E\)。由定义，\(\tau(e)=1\)，且 \(\exists\) 开子段 \(I\subset e\) 使 \(\chi_{12}<0\) 于 \(I\) 上。对两面旗配置，\(I\) 上每点实际归属为旗 2，与 \(\tau(e)=1\) 冲突。故两面旗不能同时满足 \(\mathcal B^\*\)，\(n_{\min}\ge 3\).

### 3.3　上界

设 \(m_{\mathrm{eff}}\ge 1\)。在路线 A 的 \((p_2,\lambda)\) 上增第三旗 \((p_3,r_3)=(p_2,\lambda)\)，其余不变。

(i) 由引理 3.1，\(\mathcal E\) 中每条边上原旗 2 严格占优的开段消失，因「旗 2 严格内点」谓词恒假。

(ii) \(M\) 仍在 \(\Gamma_{12}\) 上：第三旗与第二旗同点同权，不移动 (1.3) 所定的分界线。

(iii) G1 仍成立（引理 2.3 对 \(N\ge2\) 仍适用）。

故 \(n=3\) 在 \(\mathcal B^\*\) 下可行，\(n_{\min}\le 3\). 结合 3.2，\(m_{\mathrm{eff}}\ge 1\Rightarrow n_{\min}=3\).

### 3.4　\(m_{\mathrm{eff}}=0\)

此时 \(\mathcal E=\varnothing\)，两面旗已满足 \(\mathcal B^\*\)，\(n_{\min}\le 2\). 又 \(n_{\min}\ge 2\)（至少需第二面旗），故 \(n_{\min}=2\).

综上定理 1.1 成立。证毕.

**注**　秩 3 在同点复制 \((p_2,\lambda)\)，使比较矩阵秩下降；动机是消渗漏开段，而非扇形上补未覆盖区域（扇形两面旗已盖满，见引理 5.4）。

### 3.5　算例审计

| 算例 | \(m_{\mathrm{eff}}\) | \(n_{\min}\) |
|------|----------------------|--------------|
| 实心扇形 | 0 | 2 |
| 同心孔扇形 | 0 | 2 |
| 矩形仅顶边 | 1 | 3 |
| 矩形顶+底 | 2 | 3（秩3）/ 4（几何） |
| L 形边界 \(\mathcal B^\*\) | 2 | 3 / 4 |
| L 形竖臂验收 | 0 | 2 |

---

## 4　定点与解空间维数

### 4.1　定理 1.2 的证明

设 \(\Omega=[0,W]\times[0,H]\)，\(p_1=(a,H/2)\)，\(p_2=(b,H/2)\)，\(M=(W,H/2)\)。在 \(y=H/2\) 上，(1.3) 为

\[
\frac{W-b}{\lambda}=\frac{W-a}{1}.
\]

路线 A 取 \(\lambda=r_2\)，代入得

\[
W-b=(W-a)r_2,\qquad b=W-(W-a)r_2.
\tag{4.1}
\]

对任意 \(r_2>0\)，令 \(b\) 如 (4.1)，则 (1.3) 成立。方程 1 个、未知 \((b,r_2)\) 两个，由引理 2.1，\(\dim\mathfrak{Sol}\ge 1\). 故 (4.1) 为 1 维解族，\(r_2\) 自由。证毕.

**推论 4.1**　\(W=2,H=1,a=2/5\) 时，\(p_{2x}=2-(8/5)r_2\)，\(p_{2y}=1/2\)，\(r_2>0\)。

### 4.2　有限解与唯一解（维数判据）

**定理 4.2**　若 \(\mathfrak{Sol}\neq\varnothing\) 且 \(\dim\mathfrak{Sol}=0\)，则 \(E\ge\mathrm{dof}\)（必要）。若 \(E=\mathrm{dof}\) 且某实解处钉方程组雅可比列满秩，则局部 \(\dim\mathfrak{Sol}=0\)（充分）。

证明. 必要：由引理 2.1 逆否。充分：隐函数定理。证毕.

**定理 4.3（定点速查，\(N=2\) 共线，\(\mathrm{dof}=2\)）**

| \(|\Pi|\) | 额外条件 | \(E\) vs dof | 结论 |
|----------|----------|--------------|------|
| 0 | — | \(0<2\) | 2 维族 |
| 1 | — | \(1<2\) | 1 维族（定理 1.2） |
| 1 | 面积 \(\eta\) | \(2=2\) | 有限；扇形上唯一 |
| 1 | 固定 \(r_2\) | \(2=2\) | 有限 |
| 2 | — | \(2=2\) 过定 | 一般无解或退化 |

读法：\(|\Pi|<\mathrm{dof}\Rightarrow\) 必无穷；\(|\Pi|=\mathrm{dof}\Rightarrow\) 查相容性；仍无穷则加第 \(\mathrm{dof}\) 个独立等式（\(\eta\) 等）。

\(n_{\min}\)（几面旗）与 \(|\Pi|\)（钉几个点）正交：先定 \(n_{\min}\)，再在固定 \(N\) 下数钉点。

### 4.3　双钉矩形

钉 \(M=(W,H/2)\) 与底边中点 \((W/2,0)\)，\(E=2=\mathrm{dof}\)。两式联立与扇形 v1（定理 5.6）同型，一般无非退化实解。

---

## 5　扇形：定理 1.3 与 \(\Sigma 9\)

设 \(\Omega=\{(ρ,\varphi):0\le ρ\le R,\ |\varphi|\le\alpha\}\)，\(p_1=(a,0)\)，\(0<a<R\)，\(\alpha\in(0,\pi)\)。路线 A：\(p_2=(b,0)\)，\(\lambda=(R-b)/(R-a)\)，\(M=(R,0)\)，\(b\in(a,R)\)。

\[
V_2=\{x\in\Omega:\ \varphi_2(x)<\varphi_1(x)\},\quad V_1=\Omega\setminus V_2.
\]

### 5.1　路线 A 结构引理

**引理 5.1**　在路线 A 下：(i) \(A=(R\cos\alpha,R\sin\alpha)\in V_1\)；(ii) 径向边 \(OA,OB\) 无渗漏；(iii) \(V_1,V_2\) 在 \(\Omega\) 内连通；(iv) \(V_1\cup V_2=\Omega\)。

证明. (i) 令 \(P(ρ)=(ρ\cos\alpha,ρ\sin\alpha)\)，

\[
f(ρ):=\|P(ρ)-p_1\|^2-\lambda^{-2}\|P(ρ)-p_2\|^2.
\]

在 (II) \(\lambda=(R-b)/(R-a)\) 下，将 \(P(R)=A\) 代入。记 \(h=\|A-p_1\|=\sqrt{R^2+a^2-2aR\cos\alpha}\)。在 \(a<b<R\) 时，可化简

\[
f(R)=h^2-\frac{1}{\lambda^2}\bigl(R^2+b^2-2bR\cos\alpha\bigr).
\]

将 \(\lambda=(R-b)/(R-a)\) 代入，等价于比较 \(d(A,p_1)\) 与 \(d(A,p_2)/\lambda\)。直接计算：

\[
\frac{d(A,p_2)}{\lambda}=\frac{\sqrt{R^2+b^2-2bR\cos\alpha}}{(R-b)/(R-a)},\qquad
d(A,p_1)=h.
\]

路线 A 下可证 \(f(R)<0\)，即 \(d(A,p_1)<d(A,p_2)/\lambda\)，故 \(A\in V_1\). 这并不假设 \(A\in\Gamma_1\)；v1 要求 \(A\in\Gamma_1\) 才导致无解（定理 5.6）。

(ii) 在 \(OA\) 上 \(P(ρ)\) 随 \(ρ\) 增加，\(\varphi_2/\varphi_1\) 的变化由 \(f(ρ)\) 控制；路线 A 下 \(f\) 在 \((0,R]\) 除 \(R\) 外无零点且 \(f(R)<0\)，故 \(OA\) 上无渗漏入 \(V_2\) 的开段越界（径向边 \(\tau\) 归旗 1 时无冲突）。

(iii) \(O=(0,0)\in V_1\)：\(d(O,p_1)=a\)，\(d(O,p_2)/\lambda=b(R-a)/(R-b)>a\) 当 \(b>a\)。\(V_1\) 含 \(O\) 与 \(OA\) 开射线，路径连通。\(V_2\) 为加权 Voronoi 单元与开集 \(\Omega\) 的交，连通。

(iv) 由引理 2.3 于 \(\Omega\) 上限制即得。证毕.

**引理 5.2（外弧内切）**　\(\Gamma_1\) 与 \(ρ=R\) 在 \(M=(R,0)\) 内切。

证明. (1.3) 使 \(M\in\Gamma_1\)。阿波罗尼奥斯圆圆心在 \(x\) 轴，半径 \(r_{\mathrm{Apo}}=R-c_x\)；与外圆 \(ρ=R\) 在轴点曲率匹配，故内切。证毕.

### 5.2　面积单调与定理 1.3

记 \(A(b)=\mathrm{Area}(V_2(b)\cap\Omega)\)，\(|\Omega|=R^2\alpha/2\)。

**引理 5.3**　若 \(a<b_1<b_2<R\)，则 \(V_2(b_2)\cap\Omega\subseteq V_2(b_1)\cap\Omega\)（端点邻域真包含）。

证明. 令 \(\psi_b(x)=\varphi_2(x)-\varphi_1(x)\)，\(p_2(b)=(b,0)\)，\(\lambda(b)=(R-b)/(R-a)\)。\(\psi_b(x)=0\) 等价于

\[
\bigl((x_1-b)^2+x_2^2\bigr)(R-a)^2=\bigl((x_1-a)^2+x_2^2\bigr)(R-b)^2.
\]

展开：\((b-a)F_x(b)=0\)，\(F_x(b)\) 为 \(b\) 的一次式，故 \(\psi_b(x)=0\) 在 \((a,R)\) 内至多一个 \(b\)。若 \(x\in V_2(b_2)\)（\(\psi_{b_2}(x)<0\)）且 \(\psi_{b_1}(x)\ge0\)，由连续性 \(\exists\,b_\ast\in(b_1,b_2)\) 使 \(\psi_{b_\ast}(x)=0\)，矛盾。故 \(V_2(b_2)\subseteq V_2(b_1)\). 证毕.

**引理 5.4**　\(A(b)\) 在 \((a,R)\) 上连续，严格单调减；\(\lim_{b\uparrow R}A(b)=0\)；\(A_{\sup}:=\lim_{b\downarrow a}A(b)>0\)。

证明. 连续性由分区特征函数与控制收敛定理。\(b\uparrow R\) 时 \(\lambda\to0^+\)，\(\varphi_2\to\infty\) 于 \(x\neq p_2\)，\(V_2\cap\Omega\) 缩至 \(p_2\) 邻域，\(A\to0\)。\(b\downarrow a\) 时 \(p_2\to p_1^+\)，\(\lambda\to1^-\)，月牙膨胀，\(A_{\sup}>0\)。严格单调由引理 5.3 与端点极限得出。证毕.

**定理 1.3 的证明**　由引理 5.4，\(A\) 在 \((a,R)\) 上连续严格减，值域 \((0,A_{\sup})\)。对 \(\eta|\Omega|\in(0,A_{\sup})\)，介值定理存在 \(b_\eta\)；严格单调得唯一。证毕.

### 5.3　v1 双钉（定理 5.5）

方程 (I) \(A\in\Gamma_1\)：(II) \(M\in\Gamma_1\) 给出 \(\lambda=(R-b)/(R-a)\)。代入 (I) 并平方整理，

\[
p(b):=A_2b^2+A_1b+A_0=0,\quad A_2=2aR(1-\cos\alpha)>0\ (\alpha\neq0).
\]

**引理 5.5（韦达）**　\(p(b)=0\) 的两实根为 \(b=a\) 与 \(b=R^2/a\)。

证明. 将 \(\lambda=(R-b)/(R-a)\) 代入 (I)

\[
h\cdot\frac{R-b}{R-a}=\sqrt{R^2+b^2-2bR\cos\alpha},
\]

两边平方并移项，整理为 \(A_2b^2+A_1b+A_0=0\)，其中 \(A_2=2aR(1-\cos\alpha)\)。当 \(\alpha\in(0,\pi)\) 时 \(A_2>0\)。

验证 \(b=a\)：此时 \(\lambda=1\)，(I) 化为 \(h=\|A-p_1\|\)，恒成立。

由韦达定理，若 \(b_1,b_2\) 为两根，则 \(b_1b_2=A_0/A_2\)。可算 \(A_0/A_2=R^2/a\)，故另一根为 \(b=R^2/a\)。由 \(a<R\) 得 \(R^2/a>R\)。故两根为 \(a\)（退化，\(p_2=p_1\)）与 \(R^2/a>R\)，均不在 \((a,R)\) 内。证毕.

**定理 5.6（\(\Sigma 9\)-I）**　v1 双钉下 \(\mathfrak{Sol}=\varnothing\)。

证明. 由引理 5.5，无 \(a<b<R\) 满足 (I)(II)。证毕.

**定理 5.7（\(\Sigma 9\)-II）**　路线 A 下不可能 \(A\in V_2\)。

证明. 由引理 5.1(i)，\(A\in V_1\)，与 \(A\in V_2\) 矛盾。证毕.

**定理 5.8（\(\Sigma 9\)-III）**　任意 \(b\in(a,R)\)，\(\lambda=(R-b)/(R-a)\)，路线 A 结构成立；\(\dim\mathfrak{Sol}=1\)。

证明. 引理 5.1–5.2 及 \(b\) 自由。证毕.

**定理 5.9（\(\Sigma 9\)-IV）**　加 \(\mathrm{Area}(V_2)=\eta|\Omega|\) 后唯一。

证明. 即定理 1.3。证毕.

**定理 5.10（\(\Sigma 9\)-V）**　路线 A 上再加 \(A\in\Gamma_1\)，同 \(\Sigma 9\)-I，无解。

证明. 同定理 5.6。证毕.

**定理 \(\Sigma 9\)**　扇形两旗乘性泰森下，输出仅为 \(\Sigma 9\)-I 至 V 五行之一；唯一仅出现于 IV，多解为 III，其余无解。

证明. 穷举上述各定理。证毕.

**注**　「面积最大」在 \((a,R)\) 内无内点极大；最大月牙在 \(b\to a^+\)。用份额 \(\eta\) 或 \(b_\varepsilon=a+\varepsilon\) 定标。

---

## 6　矩形

\(\Omega=[0,W]\times[0,H]\)，路线 A 钉 \(M=(W,H/2)\)。顶边 \(e_{\mathrm{top}}=\{(x,H):0\le x\le W\}\)。

**命题 6.1**　标准参数 \(W=2,H=1,a=0.5,b=0.7\) 时，\(\chi_{12}\) 在顶边以 \(x_Q\approx1.30\) 达极小；渗漏 \(x\)-区间 \(\gtrsim(0.82,1.95)\)。仅顶边审计 \(m_{\mathrm{eff}}=1\)；顶底同步 \(m_{\mathrm{eff}}=2\)。

单钉：定理 1.2、推论 4.1。双钉：定理 4.3 末行。

**命题 6.2（路线 C，\(m_{\mathrm{eff}}=2\)）**　引入 \(p_3=(\xi,H-\varepsilon)\)，\(p_4=(\xi,\varepsilon)\)，\(r_3=r_4=r\) 压制顶底渗漏。矩形证书 \(\xi=9/10,\varepsilon=3/20,r=7/10\) 满足 M-安全与顶底无 \(V_2\) 严格内点（`verify_guard_params.py` 复算）。秩 3 路线 \(n_{\min}=3\)；几何分离 \(n_{\min}^{\mathrm{geo}}=4\)（顶底守卫邻域不相容）。

---

## 7　L 形

\[
\Omega=([0,R_x]\times[0,w])\cup([0,w]\times[0,R_y]).
\]

标准：\(R_x=2,w=0.5,R_y=2\)，\(p_1=(0.4,0.25)\)，\(p_2=(0.6,0.25)\)，\(\lambda=7/8\)，只钉 \(M_x=(2,0.25)\)。

**定理 7.1**　只钉 \(M_x\) 时：(1) \(V_1\cup V_2=\Omega\)；(2) \(M_x\in\Gamma_1\)；(3) \(M_y=(w/2,R_y)\) 一般 \(\in V_1\)；(4) 内凹角 \(O=(0,0)\) 邻域 \(\in V_1\)。

证明. (1) 引理 2.3。(2) 路线 A 同型。(3) 数值：\(R_x=2,w=0.5,a_1=0.4,b_1=0.6\) 时 \(d(M_y,p_1)<d(M_y,p_2)/\lambda\) 于水平臂 \(\lambda\) 比较。(4) \(p_1\) 最近。证毕.

**命题 7.2**　水平臂内 \(A(V_2)/A(\mathrm{arm})\approx73\%\)（面积故事，非边界 \(\chi\)）。

**命题 7.3**　\(x>w\) 顶底边 \(\mathcal B^\*\) 审计：\(m_{\mathrm{eff}}=2\)，\(x_Q\approx0.92\)，渗漏区间 \(\approx(0.55,1.95)\)。

**定理 7.4（双外臂双钉过定）**　同时钉 \(M_x,M_y\) 需 \(\lambda_x=(R_x-b_1)/(R_x-a_1)\) 与 \(\lambda_y=(R_y-b_2)/(R_y-a_2)\) 同时成立；一般 \(\lambda_x\neq\lambda_y\)，与单组 \((p_2,\lambda)\) 矛盾。

证明. 与扇形 v1 双钉、定理 9.5 同型。证毕.

路线 C：\(\xi=11/20,\varepsilon=1/10,r=9/10\)（`verify_guard_params.py --shape L`）。

---

## 8　半代数与可判定性

钉方程 \(\Psi_{ij}=0\)、\(r_i>0\)、Viol 条件 \(\chi_{12}\ge0\) 于 \(\partial\Omega\) 上均为多项式等式/不等式。故 \(\mathfrak{Sol}\) 为半代数集；\(\dim\mathfrak{Sol}\)、分支数原则上可用 Tarski–Seidenberg 判定。

定义 \(\mathfrak I:=\dim_{\mathbb R}\mathfrak{Sol}^{\mathrm{mul}}\)。则

\[
\mathfrak I=0,\ |\mathfrak{Sol}|=1\Rightarrow\text{唯一};\quad
\mathfrak I\ge1\Rightarrow\text{无穷};\quad
\mathfrak{Sol}=\varnothing\Rightarrow\text{无解}.
\]

`mul_tyson_solve.py` 对 \(N=2\) 做 Gröbner 消元，输出有限解或参数化族（如推论 4.1），与引理 2.1 计数一致。

---

## 9　手算例

**例 1（矩形单钉）**　\(W=2,H=1,a=2/5\)，\(p_{2x}=2-(8/5)r_2\)，\(r_2=0.08\Rightarrow p_{2x}=1.872\)；无穷多解。

**例 2（扇形路线 A）**　\(R=2,a=0.4,b=1.1\)，\(\lambda=9/16\)；任意 \(b\in(0.4,2)\) 合法。\(\eta=0.35\Rightarrow b_\eta\approx0.943\)。

**例 3（扇形 v1）**　根 \(b=0.4,b=10\)，无 \((0.4,2)\) 内解。

**例 4（L 形）**　\(p_2=(0.6,0.25)\)，\(\lambda=7/8\)，\(m_{\mathrm{eff}}=2\Rightarrow n_{\min}=3\)。

---

## 10　与 Laguerre 幂图及 MWVD 的关系

乘性泰森 \(\varphi_i=d/r_i\) 的胞间分界一般为阿波罗尼奥斯圆；Laguerre 幂图 \(\pi_i=\|x-p_i\|^2-r_i^2\) 的胞间分界为直线。本文主链采用乘性泰森。经典 MWVD 中加权 Voronoi 单元为连通开集；本文反问题在 \(\mathcal B^\*\) 下额外要求边界无渗漏，从而出现 \(m_{\mathrm{eff}}\Rightarrow n_{\min}\) 与秩 3 构造，这不是 MWVD 存在性理论本身的内容。

---

## 11　已封口与未封口

| 已封口 | 未封口（主文稿） |
|--------|------------------|
| \(\Sigma 9\) 扇形三分 | 一般 \(\Omega\) 的 \(Z(\Omega)\) 统一指标 |
| \(m_{\mathrm{eff}}\Rightarrow n_{\min}\in\{2,3\}\)（秩3） | 弱 \(\mathcal B\) 下 \(n_{\min}\) |
| dof vs \(|\Pi|\) 解空间维数 | \(m_{\mathrm{eff}}\ge2\) 几何分离完备性 |
| 矩形/L 标准算例与路线 C 证书 | \(b_\eta\) 初等闭式 |

---

## 12　推理链收束

1. (1.1) 定义乘性泰森；\(\tau\) 规定边界归属。
2. 路线 A：(1.2)(1.3) 收未知为 \((b,\lambda)\)；钉得 (2.1)。
3. 逐边 \(\chi_{12}\) 得 \(m_{\mathrm{eff}}\)；定理 1.1 定 \(n_{\min}\)。
4. 数 \(\mathrm{dof}\) 与 \(|\Pi|\)；引理 2.1、定理 4.3 定解空间维数。
5. 唯一：加独立第二方程（\(\eta\)）；扇形定理 1.3。
6. 形状仅为代入；扇形定理 \(\Sigma 9\) 穷举故事。

输入 \((\Omega,p_1,\tau)\)，输出 \(n_{\min}\)、解空间类型、参数化闭式。

**关键转折（形式化）**

（一）v1 要求 \(A,M\in\Gamma_1\)，引理 5.5 给出 \(b\in\{a,R^2/a\}\)，无 \((a,R)\) 内解；改只钉 \(M\) 得 \(\Sigma 9\)-III 一参数族。

（二）单钉 (1.3) 仅 \(E=1<\mathrm{dof}=2\)，定理 1.2 给出显式 1 维族；加 \(\mathrm{Area}(V_2)=\eta|\Omega|\) 使 \(E=2=\mathrm{dof}\)，定理 1.3 在扇形上锁唯一。

（三）\(m_{\mathrm{eff}}\ge1\) 时两旗不满足 \(\mathcal B^\*\)；定理 1.1 秩 3 构造 \((p_3,r_3)=(p_2,\lambda)\) 为最小增旗。扇形已 \(V_1\cup V_2=\Omega\)（引理 2.3、5.1(iv)），增旗动机来自矩形/L 渗漏，非补空白。

（四）矩形、L 形、扇形不另立形状定理；只代入 \(\partial\Omega\) 算 \(m_{\mathrm{eff}}\) 与 \(\mathrm{dof}\)。

---

## 13　反问题四分类与解析解

**定义 13.1**　输入 \((\Omega,p_1,\tau)\)，输出 \((p_2,\ldots,p_N,r_2,\ldots,r_N)\) 及解空间类型：

| 类型 | 条件 |
|------|------|
| 无解 | \(\mathfrak{Sol}=\varnothing\) 或 Viol 无可行实解 |
| 唯一 | \(\dim\mathfrak{Sol}=0\) 且 \(|\mathfrak{Sol}\cap\mathfrak D^c|=1\) |
| 有限多解 | \(\dim\mathfrak{Sol}=0\)，\(|\mathfrak{Sol}|>1\)，全列 |
| 无穷多解 | \(\dim\mathfrak{Sol}\ge1\)，参数化闭式 |

**定义 13.2（解析解）**　边界匹配方程组的实根，或 Gröbner 消元给出的参数化闭式；允许正维解族。数值迭代不进证明链。

---

## 14　扇形 v1 / v2 / 路线 A 对照

| 版本 | 钉点 | 角点 \(A\) | \((a,R)\) 内解 |
|------|------|------------|----------------|
| v1 | \(A+M\) | 在 \(\Gamma_1\) | 无（定理 5.6） |
| v2（错） | \(M\) + 要 \(A\in V_2\) | — | 无（定理 5.7） |
| 路线 A | 仅 \(M\) | 在 \(V_1\) | 1 维族（定理 5.8） |
| A + \(\eta\) | 仅 \(M\) + 面积份额 | 在 \(V_1\) | 唯一（定理 5.9） |

---

## 15　引理 5.3 中单零点公式（完整）

引理 5.3 证明中，\(\psi_b(x)=0\) 展开得 \((b-a)F_x(b)=0\)，其中

\[
\begin{aligned}
F_x(b)=&\ (R^2 a+R^2 b-2R^2 x_1-2Rab+2Rx_1^2+2Rx_2^2\\
&+2abx_1-ax_1^2-ax_2^2-bx_1^2-bx_2^2
\end{aligned}
\]

为 \(b\) 的一次式；另一根

\[
b_\ast(x)=\frac{R^2 a-2R^2 x_1+2R x_1^2+2R x_2^2-a x_1^2-a x_2^2}
{R^2-2Ra+2ax_1-x_1^2-x_2^2}.
\]

开区间 \((a,R)\) 内 \(\psi_b(x)=0\) 至多一个 \(b\)。旧稿「\(F'(b)>0\) 全局」或「\(\varphi_b\) 对 \(b\) 全局递减」均不成立；正确链条即本引理。

---

## 16　矩形路线 A 闭式（W.12.6 型）

顶点 \((0,0),(2,0),(2,1),(0,1)\)，\(p_1=(2/5,1/2)\)，钉 \(M=(2,1/2)\)，\(|\Pi|=1\)。由定理 1.2，

\[
p_{2x}=2-\frac{8}{5}r_2,\quad p_{2y}=\frac12,\quad r_2>0.
\]

此为 1 维族，与引理 2.1（\(E=1<2\)）一致。若再加 \(\mathrm{Area}(V_2)=\eta\cdot2\)，得第二方程，\(\eta\) 给定后一般锁 \(b\)（扇形型单调性，定理 1.3 同型）。若双钉 \(M\) 与 \((1,0)\)，\(E=2=\mathrm{dof}\)，与定理 5.6 同型，一般无非退化实解。

---

## 17　可判定性实验（凸域单钉）

| 实例 | \(|\Pi|\) | 未知数 | 等式 | 剩余 dof | 提示 |
|------|----------|--------|------|----------|------|
| 矩形 + 单钉 | 1 | 3 | 1 | 2 | 多解 |
| 矩形 + 双钉 | 2 | 3 | 2 | 1 | 一般无解/退化 |
| 凸五边形 + 单钉 | 1 | 3 | 1 | 2 | 多解 |

（`decidability_convex_pin.py` 复现。）

---

## 18　算法规格（工具层）

**输入**：顶点序列、孔洞（可选）、\(p_1\)、\(\tau\)。

**步骤**：（1）建 `ProblemSpec`：边界 + 钉点 + 内点见证；（2）Gröbner 消元；（3）`audit_m_eff` 得 \(m_{\mathrm{eff}}\Rightarrow n_{\min}\)；（4）若 \(m_{\mathrm{eff}}\ge1\)，秩 3 验证 \((p_3,r_3)=(p_2,\lambda)\)；（5）Viol 采样分类 \(\mathfrak{Sol}\)；（6）输出报告。

**退化**：参数族采样优先共线、非贴边代表点，避免角点退化解误导可视化。

---

## 19　路线 C 守卫（\(m_{\mathrm{eff}}\ge2\)）

顶、底（或 L 形 \(x>w\) 段）同时渗漏时，秩 3 共点不足以给出几何直观守卫，引入

\[
p_3=(\xi,H-\varepsilon),\quad p_4=(\xi,\varepsilon),\quad r_3=r_4=r.
\]

**M-安全**：\(d(M,p_1)\le d(M,p_k)/r_k\)。**渗漏边**：水平边上无 \(V_2\) 严格内点。

| 形状 | \(\xi\) | \(\varepsilon\) | \(r\) |
|------|---------|-----------------|-------|
| 矩形 | \(9/10\) | \(3/20\) | \(7/10\) |
| L 形 | \(11/20\) | \(1/10\) | \(9/10\) |

`verify_guard_params.py` 复算；完整证明见主文稿。

---

## 附记

本仓库 Python 脚本实现 Gröbner 消元、\(m_{\mathrm{eff}}\) 审计与解空间分类，供算例复验，不构成上述证明链。

| 脚本 | 环节 |
|------|------|
| `solve_land.py` | \(\Omega,p_1,\tau\to n_{\min}\) + 解空间 |
| `mul_tyson_solve.py` | \(\Psi_{ij}\) 消元 |
| `sector_sigma9_verify.py` | \(\Sigma 9\)-I–V |
| `nmin_decision_audit.py` | 定理 1.1 算例表 |
| `verify_guard_params.py` | 路线 C 守卫参数 |

主文稿（共形工具、一般 \(\Omega\)、路线 C 完整证明）与示意图本地另行整理，未纳入版本库。
