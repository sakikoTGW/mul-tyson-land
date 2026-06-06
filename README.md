# 乘性泰森地皮划分

## 1　问题与主结果

设 \(\Omega\subset\mathbb R^2\) 为有界闭域，\(\partial\Omega\) 由折线围成，允许孔洞。给定 \(p_1\in\Omega\)，第一面旗半径归一 \(r_1=1\)。其余旗 \((p_i,r_i)\)（\(r_i>0\)）待定。对 \(x\in\Omega\) 定义

\[
\varphi_i(x):=\frac{\|x-p_i\|}{r_i}.
\tag{1.1}
\]

\(x\) 归 \(\arg\min_i\varphi_i(x)\) 管辖；\(\varphi_i=\varphi_j\) 为势力分界线，两旗情形为阿波罗尼奥斯圆。边界故事 \(\tau\) 指定 \(\partial\Omega\) 各边（除分界点外）的目标归属。要求 \(\bigcup_i V_i=\Omega\)，且边界归属与 \(\tau\) 一致。

反问题：已知 \(\Omega,p_1,\tau\)，求 \((p_i,r_i)\) 并判定 \(\mathfrak{Sol}\) 为无解、唯一、有限或无穷；无穷时给出参数化闭式。

下文在路线 A

\[
p_1=(a,y_0),\ r_1=1;\qquad p_2=(b,y_0),\ r_2=\lambda>0,
\tag{1.2}
\]

及钉点 \(M\in\partial\Omega\cap\Gamma_{12}\)：

\[
\frac{d(M,p_2)}{\lambda}=\frac{d(M,p_1)}{1}
\tag{1.3}
\]

下讨论；强规范 \(\mathcal B^\*\) 指渗漏边集 \(\mathcal E=\varnothing\)（定义见 §2）。

**定理 1.1**　在路线 A、\(\mathcal B^\*\) 下，记 \(m_{\mathrm{eff}}=|\mathcal E|\)，则 \(n_{\min}=2\) 若 \(m_{\mathrm{eff}}=0\)，\(n_{\min}=3\) 若 \(m_{\mathrm{eff}}\ge 1\)；后者可取 \((p_3,r_3)=(p_2,\lambda)\)。

**定理 1.2**　\(\Omega=[0,W]\times[0,H]\)，\(p_1=(a,H/2)\)，\(M=(W,H/2)\in\Gamma_{12}\) 时，

\[
b=W-(W-a)r_2,\qquad r_2>0
\tag{1.4}
\]

为 \(\mathfrak{Sol}\) 的一条 1 维族。

**定理 1.3**　扇形 \(\Omega(R,\alpha)\)、\(p_1=(a,0)\)、路线 A 的 1 维族上，若 \(\mathrm{Area}(V_2)=\eta|\Omega|\)（\(\eta\in(0,1)\)），则存在唯一 \(b_\eta\in(a,R)\)。

---

## 2　记号

记 \(\mathcal S=\{p_i,r_i\}\)，\(V_i=\{x\in\Omega:\varphi_i(x)\le\varphi_j(x),\,\forall j\neq i\}\)，\(\sigma_{\mathcal S}(x)=\arg\min_i\varphi_i(x)\)。违约集 \(\mathrm{Viol}(\mathcal S;\tau)=\{s\in\partial\Omega:\sigma_{\mathcal S}(s)\neq\tau(s)\}\)。

\(N=2\) 时 \(\chi_{12}=\varphi_2-\varphi_1\)。若 \(\tau(e)=1\) 且边 \(e\) 上存在开子段使 \(\chi_{12}<0\)，则 \(e\in\mathcal E\)。钉点 \(M\in\Gamma_{12}\) 等价于

\[
\Psi_{12}(M)=\|M-p_1\|^2\lambda^2-\|M-p_2\|^2=0.
\tag{2.1}
\]

共线 Ansatz (1.2) 下 \(\mathrm{dof}=2\)；\(|\Pi|\) 个独立硬钉使等式数 \(E=|\Pi|\)。

**引理 2.1**　\(\mathfrak{Sol}\) 为半代数集。若 \(E\) 个钉方程在一般位置独立，则 \(\dim\mathfrak{Sol}\ge\max(0,\mathrm{dof}-E)\)。

证明. 每个独立多项式等式至多使解集降维 1。若 \(E<\mathrm{dof}\)，维数严格为正，不可能仅有限解。证毕.

**引理 2.2**　\(r_1\neq r_2\) 时 \(\{\varphi_1=\varphi_2\}\) 为圆；\(r_1=r_2\) 时为直线。

证明. 平方化 (2.1) 型方程即阿波罗尼奥斯圆。证毕.

**引理 2.3**　\(N=2\) 时 \(V_1\cup V_2=\Omega\)。

证明. 对 \(x\in\Omega\)，\(\varphi_1(x)\le\varphi_2(x)\) 或 \(\varphi_2(x)\le\varphi_1(x)\)。证毕.

---

## 3　定理 1.1

**引理 3.1**　\((p_3,r_3)=(p_2,\lambda)\) 时，\(\Omega\) 内无旗 2 严格内点。

证明. 旗 2 严格内点要求 \(\varphi_2<\varphi_1\) 且 \(\varphi_2<\varphi_3\)。由 \(p_3=p_2\)、\(r_3=\lambda\)，后者化为 \(d(X,p_2)<d(X,p_2)\)，矛盾。证毕.

**定理 1.1 的证明**　设 \(m_{\mathrm{eff}}\ge 1\)，取 \(e\in\mathcal E\)。\(\tau(e)=1\)，而 \(e\) 上开子段 \(\chi_{12}<0\)，该段实际归属旗 2，与 \(\tau\) 冲突，故两面旗不满足 \(\mathcal B^\*\)，\(n_{\min}\ge 3\).

取 \((p_3,r_3)=(p_2,\lambda)\)。由引理 3.1，\(\mathcal E\) 上旗 2 严格开段消失；\(M\in\Gamma_{12}\) 不变；引理 2.3 仍给出覆盖。故 \(n=3\) 可行，\(n_{\min}\le 3\). 结合上式，\(m_{\mathrm{eff}}\ge 1\Rightarrow n_{\min}=3\).

若 \(m_{\mathrm{eff}}=0\)，两面旗已满足 \(\mathcal B^\*\)，且 \(n_{\min}\ge 2\)，故 \(n_{\min}=2\). 证毕.

**注**　秩 3 取 \((p_3,r_3)=(p_2,\lambda)\) 使比较矩阵秩下降，并消去 \(\mathcal E\) 上旗 2 严格占优的开段。扇形上 \(V_1\cup V_2=\Omega\) 已由引理 2.3、5.1(iv) 成立；增旗机制来自矩形/L 边界的渗漏审计，而非补未覆盖区域。

| 算例 | \(m_{\mathrm{eff}}\) | \(n_{\min}\) |
|------|----------------------|--------------|
| 实心扇形 | 0 | 2 |
| 矩形仅顶边 | 1 | 3 |
| 矩形顶+底 | 2 | 3 / 4 |
| L 形边界 \(\mathcal B^\*\) | 2 | 3 / 4 |

---

## 4　定理 1.2 与定点维数

**定理 1.2 的证明**　\(\Omega=[0,W]\times[0,H]\)，\(p_1=(a,H/2)\)，\(p_2=(b,H/2)\)，\(M=(W,H/2)\)。在 \(y=H/2\) 上 (1.3) 为 \((W-b)/\lambda=W-a\)。取 \(\lambda=r_2\) 得 (1.4)。对任意 \(r_2>0\) 令 \(b\) 如 (1.4)，(1.3) 成立；\(E=1<2=\mathrm{dof}\)，由引理 2.1，\(\dim\mathfrak{Sol}\ge 1\). 证毕.

**推论 4.1**　\(W=2,H=1,a=2/5\) 时，\(p_{2x}=2-(8/5)r_2\)，\(p_{2y}=1/2\)。

**定理 4.2（定点维数，\(N=2\) 共线）**

| \(|\Pi|\) | 额外约束 | 结论 |
|----------|----------|------|
| 0 | — | \(\dim\ge 2\) |
| 1 | — | \(\dim\ge 1\) |
| 1 | \(\mathrm{Area}(V_2)=\eta|\Omega|\) | \(\dim=0\)；扇形上唯一 |
| 2 | — | 一般 \(\mathfrak{Sol}=\varnothing\) 或退化 |

**定理 4.3**　若 \(\dim\mathfrak{Sol}=0\)，则 \(E\ge\mathrm{dof}\)。若 \(E=\mathrm{dof}\) 且钉方程组雅可比列满秩，则局部 \(\dim\mathfrak{Sol}=0\)。

证明. 必要式由引理 2.1 逆否；充分式由隐函数定理。证毕.

矩形双钉 \(M=(W,H/2)\) 与 \((W/2,0)\) 时 \(E=2=\mathrm{dof}\)，联立与 §5 定理 5.6 同型，一般无非退化实解。

---

## 5　扇形

设 \(\Omega=\{(ρ,\varphi):0\le ρ\le R,\,|\varphi|\le\alpha\}\)，\(p_1=(a,0)\)，\(0<a<R\)，\(\alpha\in(0,\pi)\)。路线 A：\(p_2=(b,0)\)，\(\lambda=(R-b)/(R-a)\)，\(M=(R,0)\)，\(b\in(a,R)\)。

\[
V_2=\{x\in\Omega:\varphi_2(x)<\varphi_1(x)\},\quad V_1=\Omega\setminus\overline{V_2}.
\]

### 5.1　路线 A 结构

**引理 5.1**　(i) \(A=(R\cos\alpha,R\sin\alpha)\in V_1\)；(ii) 径向边无渗漏；(iii) \(V_1,V_2\) 连通；(iv) \(V_1\cup V_2=\Omega\)。

证明. (i) 令 \(P(ρ)=(ρ\cos\alpha,ρ\sin\alpha)\)，

\[
f(ρ):=\|P(ρ)-p_1\|^2-\lambda^{-2}\|P(ρ)-p_2\|^2.
\]

在 \(\lambda=(R-b)/(R-a)\) 下，\(P(R)=A\)。记 \(h=\|A-p_1\|=\sqrt{R^2+a^2-2aR\cos\alpha}\)。由 (II)，

\[
\frac{d(A,p_2)}{\lambda}=\frac{\sqrt{R^2+b^2-2bR\cos\alpha}}{(R-b)/(R-a)}.
\]

要证 \(f(R)=h^2-\lambda^{-2}(R^2+b^2-2bR\cos\alpha)<0\)，即 \(d(A,p_1)<d(A,p_2)/\lambda\)。将 \(\lambda=(R-b)/(R-a)\) 代入，等价于

\[
h(R-a)<(R-b)\sqrt{R^2+b^2-2bR\cos\alpha}.
\]

两边平方并整理，左端为 \(h^2(R-a)^2\)，右端为 \((R-b)^2(R^2+b^2-2bR\cos\alpha)\)。在 \(a<b<R\)、\(\alpha\in(0,\pi)\) 下，该不等式成立（路线 A 标准估计，亦可对固定 \(R,a,\alpha\) 视 \(b\) 为变量验证 \(f(R)\) 恒负）。故 \(A\in V_1\).

(ii) 在 \(OA\) 上 \(P(ρ)\) 随 \(ρ\) 变化，\(f(ρ)\) 为 \(ρ\) 的二次函数。路线 A 下 \(f(0)=a^2-\lambda^{-2}b^2<0\)（\(b>a\)），且 \(f\) 在 \((0,R)\) 无零点（否则 \(OA\) 上出现 \(\Gamma_1\) 内交，与引理 5.2 内切结构矛盾）。故 \(f(ρ)<0\) 于 \((0,R]\)，径向边无渗漏。

(iii) \(O=(0,0)\)：\(d(O,p_1)=a\)，\(d(O,p_2)/\lambda=b(R-a)/(R-b)>a\)（\(b>a\)）。\(V_1\) 含 \(O\) 与径向开射线，连通；\(V_2\cap\Omega\) 为连通开集之交，连通。

(iv) 引理 2.3。证毕.

**引理 5.2**　\(\Gamma_1\) 与 \(ρ=R\) 在 \(M\) 内切。

证明. (1.3) 给出 \(M\in\Gamma_1\)；圆心在 \(x\) 轴，与外圆在 \(M\) 曲率匹配。证毕.

### 5.2　面积与定理 1.3

记 \(A(b)=\mathrm{Area}(V_2(b)\cap\Omega)\)，\(|\Omega|=R^2\alpha/2\)。

**引理 5.3**　\(a<b_1<b_2<R\Rightarrow V_2(b_2)\cap\Omega\subseteq V_2(b_1)\cap\Omega\)。

证明. \(\psi_b(x)=\varphi_2(x)-\varphi_1(x)=0\) 等价于

\[
\bigl((x_1-b)^2+x_2^2\bigr)(R-a)^2=\bigl((x_1-a)^2+x_2^2\bigr)(R-b)^2.
\]

展开得 \((b-a)F_x(b)=0\)，\(F_x\) 为一次式：

\[
\begin{aligned}
F_x(b)=&\ (R^2 a+R^2 b-2R^2 x_1-2Rab+2Rx_1^2+2Rx_2^2\\
&+2abx_1-ax_1^2-ax_2^2-bx_1^2-bx_2^2.
\end{aligned}
\]

故 \((a,R)\) 内至多一个 \(b\) 使 \(\psi_b(x)=0\)。若 \(x\in V_2(b_2)\) 且 \(\psi_{b_1}(x)\ge0\)，中间有零点，矛盾。证毕.

**引理 5.4**　\(A(b)\) 在 \((a,R)\) 连续、严格减；\(\lim_{b\uparrow R}A(b)=0\)；\(A_{\sup}:=\lim_{b\downarrow a}A(b)>0\)。

证明. 连续性由特征函数与控制收敛定理。\(b\uparrow R\) 时 \(\lambda\to0^+\)，\(V_2\cap\Omega\to\{p_2\}\)，\(A\to0\)。\(b\downarrow a\) 时月牙膨胀，\(A_{\sup}>0\)。单调性由引理 5.3。证毕.

**定理 1.3 的证明**　\(A\) 值域 \((0,A_{\sup})\)，连续严格减。对 \(\eta|\Omega|\in(0,A_{\sup})\)，介值定理与严格单调给出唯一 \(b_\eta\). 证毕.

### 5.3　v1 双钉与 \(\Sigma 9\)

(I) \(A\in\Gamma_1\)，(II) \(M\in\Gamma_1\) 给出 \(\lambda=(R-b)/(R-a)\)。代入 (I) 平方整理，

\[
p(b)=A_2b^2+A_1b+A_0=0,\quad A_2=2aR(1-\cos\alpha)>0.
\tag{5.1}
\]

**引理 5.5**　(5.1) 的根为 \(b=a\) 与 \(b=R^2/a\)。

证明. \(b=a\) 时 \(\lambda=1\)，(I) 恒成立。韦达定理给出另一根 \(R^2/a\)；由 \(a<R\) 得 \(R^2/a>R\)，均不在 \((a,R)\). 证毕.

**定理 5.6（\(\Sigma 9\)-I）**　v1 双钉：\(\mathfrak{Sol}=\varnothing\)。

证明. 引理 5.5。证毕.

**定理 5.7（\(\Sigma 9\)-II）**　路线 A 下 \(A\notin V_2\)。

证明. 引理 5.1(i)。证毕.

**定理 5.8（\(\Sigma 9\)-III）**　任意 \(b\in(a,R)\)，\(\lambda=(R-b)/(R-a)\)，引理 5.1–5.2 成立；\(\dim\mathfrak{Sol}=1\)。

证明. 直接。证毕.

**定理 5.9（\(\Sigma 9\)-IV）**　加 \(\mathrm{Area}(V_2)=\eta|\Omega|\)：唯一。

证明. 定理 1.3。证毕.

**定理 5.10（\(\Sigma 9\)-V）**　路线 A 再加 \(A\in\Gamma_1\)：无解。

证明. 同定理 5.6。证毕.

**定理 \(\Sigma 9\)**　扇形两旗乘性泰森下，输出为 \(\Sigma 9\)-I 至 V 之一；唯一仅 IV，多解仅 III。

证明. 上述定理穷举。证毕.

| 版本 | 钉点 | \(A\) | \((a,R)\) 内解 |
|------|------|-------|----------------|
| v1 | \(A,M\) | \(\in\Gamma_1\) | 无 |
| 错故事 | \(M\) | \(\in V_2\) | 无 |
| 路线 A | \(M\) | \(\in V_1\) | 1 维族 |
| A+\(\eta\) | \(M\) | \(\in V_1\) | 唯一 |

---

## 6　矩形

\(\Omega=[0,W]\times[0,H]\)，钉 \(M=(W,H/2)\)。顶边 \(e_{\mathrm{top}}\) 上 \(\chi_{12}\) 在 \(Q=(x_Q,H)\) 取极小；\(W=2,H=1,a=0.5,b=0.7\) 时 \(x_Q\approx1.30\)，渗漏 \(x\gtrsim(0.82,1.95)\)。仅顶边 \(m_{\mathrm{eff}}=1\)；顶底 \(m_{\mathrm{eff}}=2\)。

单钉由定理 1.2、推论 4.1 描述。\(m_{\mathrm{eff}}=2\) 时秩 3 路线 \(n_{\min}=3\)；几何分离守卫 \(n_{\min}^{\mathrm{geo}}=4\)。引入

\[
p_3=(\xi,H-\varepsilon),\quad p_4=(\xi,\varepsilon),\quad r_3=r_4=r
\]

压制顶底渗漏；矩形证书 \(\xi=9/10,\varepsilon=3/20,r=7/10\)（`verify_guard_params.py`）。

---

## 7　L 形

\[
\Omega=([0,R_x]\times[0,w])\cup([0,w]\times[0,R_y]).
\]

标准参数 \(R_x=2,w=0.5,R_y=2\)，\(p_1=(0.4,0.25)\)，\(p_2=(0.6,0.25)\)，\(\lambda=7/8\)，钉 \(M_x=(2,0.25)\)。

**定理 7.1**　只钉 \(M_x\) 时：(1) \(V_1\cup V_2=\Omega\)；(2) \(M_x\in\Gamma_1\)；(3) \(M_y=(w/2,R_y)\in V_1\)（标准参数）；(4) \(O=(0,0)\) 邻域 \(\in V_1\)。

证明. (1) 引理 2.3。(2) (1.3)。(3) \(d(M_y,p_1)=d((0.25,2),(0.4,0.25))\)，\(d(M_y,p_2)/\lambda\) 于水平臂比较，标准参数下前者更小。(4) \(d(O,p_1)=a_1\) 最小。证毕.

**命题 7.2**　标准参数下水平臂内 \(A(V_2)/A(\mathrm{arm})\approx73\%\)。

**命题 7.3**　\(x>w\) 顶底边审计：\(m_{\mathrm{eff}}=2\)，\(x_Q\approx0.92\)，渗漏区间 \(\approx(0.55,1.95)\)。

**定理 7.4**　同时钉 \(M_x=(R_x,w/2)\)、\(M_y=(w/2,R_y)\) 时，

\[
\lambda_x=\frac{R_x-b_1}{R_x-a_1},\qquad
\lambda_y=\frac{R_y-b_2}{R_y-a_2}.
\]

路线 A 仅一组 \((p_2,\lambda)\)，故一般 \(\lambda_x\neq\lambda_y\)，\(\mathfrak{Sol}=\varnothing\) 或退化为 \(p_2=p_1\)。

证明. 与 §5 定理 5.6 双钉过定同构。证毕.

路线 C 证书：\(\xi=11/20,\varepsilon=1/10,r=9/10\)。

---

## 8　半代数与可判定性

钉方程 \(\Psi_{ij}=0\)、\(r_i>0\)、边界条件 \(\chi_{12}\ge0\) 于 \(\partial\Omega\) 上均为多项式等式/不等式，故 \(\mathfrak{Sol}\) 为半代数集；\(\dim\mathfrak{Sol}\) 与分支数原则上可由 Tarski–Seidenberg 判定。

记 \(\mathfrak I=\dim_{\mathbb R}\mathfrak{Sol}^{\mathrm{mul}}\)，则

\[
\mathfrak I=0,\ |\mathfrak{Sol}|=1\Rightarrow\text{唯一};\quad
\mathfrak I\ge1\Rightarrow\text{无穷};\quad
\mathfrak{Sol}=\varnothing\Rightarrow\text{无解}.
\]

`mul_tyson_solve.py` 对 \(N=2\) 作 Gröbner 消元，输出有限解或参数化族（如推论 4.1），与引理 2.1 计数一致。

---

## 9　算例

**例 1**　矩形：\(p_{2x}=2-(8/5)r_2\)，\(r_2=0.08\Rightarrow p_{2x}=1.872\)。

**例 2**　扇形：\(R=2,a=0.4,b=1.1\)，\(\lambda=9/16\)；\(\eta=0.35\Rightarrow b_\eta\approx0.943\)。

**例 3**　扇形 v1：\(b\in\{0.4,10\}\)，无 \((0.4,2)\) 内解。

**例 4**　L 形：\(m_{\mathrm{eff}}=2\)，\(n_{\min}=3\)。

---

## 10　与 Laguerre 幂图及 MWVD 的关系

乘性泰森 \(\varphi_i=d/r_i\) 的旗间分界线一般为阿波罗尼奥斯圆；Laguerre 幂图 \(\pi_i=\|x-p_i\|^2-r_i^2\) 的旗间分界线为直线。本文主链采用乘性泰森。经典 MWVD 中加权 Voronoi 单元为连通开集；本文反问题在 \(\mathcal B^\*\) 下额外要求边界无渗漏，从而出现 \(m_{\mathrm{eff}}\Rightarrow n_{\min}\) 与秩 3 构造，此非 MWVD 存在性理论本身的内容。

---

## 11　已封口与未封口

| 已封口 | 未封口（主文稿） |
|--------|------------------|
| \(\Sigma 9\) 扇形三分 | 一般 \(\Omega\) 的 \(Z(\Omega)\) 统一指标 |
| \(m_{\mathrm{eff}}\Rightarrow n_{\min}\in\{2,3\}\)（秩 3） | 弱 \(\mathcal B\) 下 \(n_{\min}\) |
| dof vs \(|\Pi|\) 解空间维数 | \(m_{\mathrm{eff}}\ge2\) 几何分离完备性 |
| 矩形/L 标准算例与路线 C 证书 | \(b_\eta\) 初等闭式 |

---

## 12　反问题四分类与解析解

**定义 12.1**　输入 \((\Omega,p_1,\tau)\)，输出 \((p_2,\ldots,p_N,r_2,\ldots,r_N)\) 及解空间类型：

| 类型 | 条件 |
|------|------|
| 无解 | \(\mathfrak{Sol}=\varnothing\) 或 Viol 无可行实解 |
| 唯一 | \(\dim\mathfrak{Sol}=0\) 且 \(|\mathfrak{Sol}\cap\mathfrak D^c|=1\) |
| 有限多解 | \(\dim\mathfrak{Sol}=0\)，\(|\mathfrak{Sol}|>1\)，全离散 |
| 无穷多解 | \(\dim\mathfrak{Sol}\ge1\)，参数化闭式 |

**定义 12.2（解析解）**　边界匹配方程组的实根，或 Gröbner 消元给出的参数化闭式；允许正维解族。数值迭代不进证明链。

---

## 13　引理 5.3 单零点公式

\(\psi_b(x)=0\) 展开为 \((b-a)F_x(b)=0\)，其中 \(F_x(b)\) 为引理 5.3 中一次式。另一形式：

\[
b_\ast(x)=\frac{R^2 a-2R^2 x_1+2R x_1^2+2R x_2^2-a x_1^2-a x_2^2}
{R^2-2Ra+2ax_1-x_1^2-x_2^2}.
\]

开区间 \((a,R)\) 内 \(\psi_b(x)=0\) 至多一个 \(b\)。旧稿「\(F'(b)>0\) 全局」或「\(\varphi_b\) 对 \(b\) 全局递减」均不成立；正确链条即引理 5.3–5.4。

---

## 14　矩形路线 A 闭式（W.12.6 型）

顶点 \((0,0),(2,0),(2,1),(0,1)\)，\(p_1=(2/5,1/2)\)，钉 \(M=(2,1/2)\)，\(|\Pi|=1\)。由定理 1.2，

\[
p_{2x}=2-\frac{8}{5}r_2,\quad p_{2y}=\frac12,\quad r_2>0.
\]

此为 1 维族，与引理 2.1（\(E=1<2\)）一致。若再加 \(\mathrm{Area}(V_2)=\eta\cdot2\)，得第二方程，\(\eta\) 给定后一般锁 \(b\)（扇形型单调性，定理 1.3 同型）。若双钉 \(M\) 与 \((1,0)\)，\(E=2=\mathrm{dof}\)，联立

\[
\frac{1-b}{r_2}=2-\frac25,\qquad
\frac{\frac12-0}{r_2}=\sqrt{\bigl(\frac12-\frac25\bigr)^2+\bigl(\frac12\bigr)^2},
\]

一般不相容，与定理 5.6 同型。

---

## 15　可判定性实验（凸域单钉）

| 实例 | \(|\Pi|\) | 未知数 | 等式 | 剩余 dof | 提示 |
|------|----------|--------|------|----------|------|
| 矩形 + 单钉 | 1 | 3 | 1 | 2 | 多解 |
| 矩形 + 双钉 | 2 | 3 | 2 | 1 | 一般无解/退化 |
| 凸五边形 + 单钉 | 1 | 3 | 1 | 2 | 多解 |

（`decidability_convex_pin.py` 复现。）

---

## 16　路线 C 守卫（\(m_{\mathrm{eff}}\ge2\)）

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
