# Project 2 研究报告

Author: Yigu Wang

Date: 2026.01.06

## 选择

我选了Cusp Potential作为研究对象，因为它在0点有一个尖锐的拐角，这种奇异结构会对波函数的行为产生显著影响。Cusp Potential的数学形式为：
$$
V(x) = V_0 |x|^{\alpha}, \quad 0 < \alpha < 1
$$

## Part A — 本征态与波函数结构

首先需要求解本征态。我用有限差分的方法，把一维定态薛定谔方程离散化成一个特征值问题。代码在 `quantum_tunneling/bound_states.py` 里，就是搭一个稀疏矩阵然后用 scipy 求特征值。

### 设置与求解过程

打开 `part_a_bound_states.ipynb`，配置好网格参数（`L` 和 `N`），选择想研究的势函数（我选的是 Cusp，但也可以试试 Exponential well 或 Soft-barrier）。运行 `run_bound_states(cfg)` 就能得到一组本征能和对应的波函数。

得到的结果里有：

- `E`：本征能量，按升序排好
- `psi`：波函数，列向量对应不同的量子态
- `metrics`：每个态的一些指标，比如均值、方差、IPR（越大越局域）
- `forbidden`：禁阻区概率，就是波函数在 V>E 地方的积分

### 看势能和波函数的样子

先画个图看看势能长什么样，以及几个本征态分别在哪个能级：

![FIGURE/part_a_plot_potential_and_states](figures/part_a_plot_potential_and_states.png)

从图上能看出 Cusp 势在 x=0 处有个尖点，波函数在低能级时比较集中，高能级时开始向外伸。

### 改变势深看会发生什么

想理解势形状对本征态的影响，我对 V0 参数做了一个扫描。从最弱的势扫到很深的势，看基态的 σ（展宽）和 IPR（局域性）怎么变化。

跑了一堆计算之后发现：势越深，基态越被"压"进中心，σ 变小，IPR 变大。这是直观的 —— 更陡峭的阱就更容易把粒子困在原地。

![FIGURE/ground_state_evolution.gif](figures/ground_state_evolution.gif)

这个动画展示了基态随 V0 的演化。看得出随着 V0 增大，波函数确实在逐渐收缩。

### 被禁阻区吸引的概率

再看看对于不同的量子态 n，禁阻区的概率是多少。这反映波函数向势垒外"泄露"的程度。

结果比较有趣：能级越高，禁阻区概率反而越小。我想这是可能因为Vcusp在远离中心时增长得很快，导致高能态的波函数主要集中在势阱内。

![FIGURE/part_a_forbidden_region_probability_vs_n](figures/part_a_forbidden_region_probability_vs_n.png)

### 波函数伸出来多远

定义一个"穿透深度"：从转折点开始，波函数的概率密度衰减到 1/e 时的距离。计算这个深度对不同 n 的变化。

我发现：随着 n 增大，穿透深度也增大。我认为这是因为随着能级升高，Vcusp开始变得平缓，波函数可以更有效地伸展出去。

![FIGURE/part_a_forbidden_penetration_depth_vs_n](figures/part_a_forbidden_penetration_depth_vs_n.png)

再来看看对于不同 n，波函数在禁阻区的空间衰减率是如何变化的。

![FIGURE/part_a_forbidden_penetration_decay_rate.png](figures/part_a_forbidden_penetration_decay_rate.png)

我同样认为这是因为随着能级升高，Vcusp平缓，波函数可以在禁阻区延伸受到的阻碍更小。

## Part C — 外场下的势垒变形与隧穿性质

在这一部分，我加上一个均匀外电场 $F$，研究 exotic barriers 在倾斜下如何变形，以及何时会发生过势垒逃逸。

### 方法

对原势施加线性倾斜：$V(x) \to V(x) - Fx$。扫描一系列外场强度 $F$，对每个 $F$ 值用 WKB 方法计算透射率。代码实现在 `part_c_field_scan.ipynb` 中，核心是调用 `run_field_scan` 进行多点计算。

### 势垒随外场的变形

首先看势能图像怎么变。当 $F$ 增大时，右侧势垒被逐渐压低，直到某个临界值时，势垒顶部下降到能级以下，粒子就可以直接逃逸。

![FIGURE: part_c_tilted_barrier_animation](figures/part_c_tilted_barrier_animation.gif)

从动画可以看出，Cusp 势在外场作用下逐渐被"拉斜"。$x=0$处的尖角依然存在（这是 singular 的体现），但势垒顶部的高度在快速下降。

### 透射率随外场的变化

现在看 WKB 透射率 $T_{WKB} = \exp(-2S/\hbar)$ 怎么随 $F$ 变化。

![FIGURE: part_c_transmission_vs_field_ln](figures/part_c_transmission_vs_field_ln.png)

在 $\ln T$ 图上，随着 $F$ 增大，$\ln T$ 逐渐上升，表明透射率在增加。注意到在某个 $F$ 值($F \approx 1.87$)之后，$\ln T = 0$，这意味着越垒电离将会发生，隧穿概率变为1。

![FIGURE: part_c_transmission_vs_field_linear](figures/part_c_transmission_vs_field_linear.png)

从线性图更清楚地看到：透射率在某个 $F$ 值之后突然增长。这个转折对应于 ultrasoft barrier 向 over-the-barrier 的过渡。最后在某个 $F$ 值($F \approx 1.87$)之后，越垒电离发生，透射率达到1。

### 过势垒逃逸的临界场强

对 Cusp 势 $V(x) = V_0|x|^{\alpha}$，施加外场后势垒消失的条件是势顶 $V_{top} = E_n$。通过求导找势顶位置，再代回求解，可以得到临界场强：

$$
F_{ionization} = \alpha \left( \frac{(1-\alpha)^{(1-\alpha)} V_0}{ E_n^{(1-\alpha)}} \right)^{1/\alpha}
$$

实现这个公式在 notebook 第 7 个代码单元。

### 不同量子态的逃逸阈值

关键发现：**高能态更容易逃逸**。

![FIGURE: part_c_ionization_field_vs_state](figures/part_c_ionization_field_vs_state.png)

从图上看，$F_{ionization}$ 随着能级 $n$ 增加而单调下降。物理直觉是：粒子能量越高，用越少的外界帮助就能越过势垒。

具体数据来自第 8 个代码单元的输出，每个激发态都对应一个不同的 $F_{ionization}$ 值。

这意味着如果我要通过外场让某个特定的激发态逃逸，需要的场强比基态小得多。我认为这在和原子电离相关的实验上有重要意义。

### Ultrasoft 和 Singular 势垒的具体条件

**Ultrasoft barrier 的形成条件：**

势垒变得 ultrasoft 当且仅当：
$$
V_{top}(F) \approx E_n, \quad \text{即} \quad F \lesssim F_{ionization}
$$

在 $F_{ionization} - F \ll F_{ionization}$ 时，势垒高度 $\Delta V = V_{top} - E_n$ 非常小。此时：

- 作用量积分 $S = \int_{x_1}^{x_2} \sqrt{2m(V(x) - E_n)} dx$ 快速趋向于零
- WKB 透射率 $T_{WKB} = \exp(-2S/\hbar)$ 快速趋向于 1
- 但这不是真正的"隧穿增强"，而是势垒本身消失了，WKB 公式失效

在 Cusp 势的情况下，ultrasoft barrier 出现在：
$$
F \approx F_{ionization} = \alpha \left( \frac{(1-\alpha)^{(1-\alpha)} V_0}{ E_n^{(1-\alpha)}} \right)^{1/\alpha}
$$

从我的数据看，基态的 $F_{ionization} \approx 1.87$，所以当 $F > 1.8$ 时势垒已经很平缓了。

**Singular barrier 的形成条件：**

Singular 不是由外场引起的，而是势函数本身的固有性质。对 Cusp 势：
$$
V(x) = V_0|x|^{\alpha}, \quad 0 < \alpha < 1
$$

在 $x = 0$ 处的导数：
$$
\frac{dV}{dx} = \alpha V_0 \text{sgn}(x) |x|^{\alpha-1} \to \infty \quad \text{当} \, x \to 0
$$

这个无穷导数（或"尖点"）即使加上外场也不会消失：
$$
V(x) - Fx = V_0|x|^{\alpha} - Fx
$$

在 $x=0$ 处仍然有尖点。因此 Singular barrier 的条件就是：**Cusp 势的指数 $\alpha < 1$ 时，尖点总是存在**。

这对数值计算有影响：

- 当转折点 $x_1$ 或 $x_2$ 离 $x=0$ 很近时，$\sqrt{V(x) - E}$ 的陡峭性会导致 WKB 积分数值不稳定
- 二阶 WKB 修正（涉及势的二阶导数）在尖点处发散

**两者的关系：**

Ultrasoft 和 Singular 是独立的两个问题：

- **Ultrasoft** 是外场导致的、可控的（通过改变 $F$ 大小），在越过 $F_{ionization}$ 后消失
- **Singular** 是势形本身的特性、不可避免的，对所有 $F$ 都存在

## Part D — 时间依赖薛定谔动力学（TDSE）

在 Part D 我关注时间动力学，直接利用Split Operator + FFT方法数值求解波函数的时间演化，从而提取逃逸率与通量信息。这里使用的一维 TDSE 为：

$$
i\hbar\,\frac{\partial \psi(x,t)}{\partial t}
=\left[-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2}+V(x)-Fx\right]\psi(x,t)
$$

对应的主流程写在 `part_d_tdse_v2.ipynb`：先求定态本征态作为初态，再把势能替换为带外场倾斜的 $V(x)-Fx$，最后用Split Operator + FFT方法演化波函数。

### 一组具体参数

我先固定一个势函数和网格参数：

- 势：Cusp，$V(x)=V_0|x|^{\alpha}$，例如 `V0=10.0, alpha=0.5`
- 空间网格：`L=100.0, N=1200`
- 外场：例如 `F=1.5`
- 时间步进：`dt=0.002, duration=64.0, record_interval=25`
- 吸收边界（CAP）：在边界附近加一个复吸收项，避免反射污染（配置在 `cfg['tdse']['cap']`）

这一块的配置和参数都在 `part_d_tdse_v2.ipynb` 的第 2 个代码单元中， `cfg = {...}` 里。CAP 的主要目的，是在波函数到达边界前将其吸收，避免边界反射回到阱内区域，从而污染逃逸率的估计。

由以上参数形成的倾斜势能和初始态波函数如下图所示：

![FIGURE: part_d_tilted_potential_and_bound_states](figures/part_d_tilted_potential_and_bound_states.png)

### 用阱内概率定义逃逸 $\Gamma_n$

PS. 这一块我只关注一个特定的本征态`state_index`=3（第4个能级）。运行 TDSE 得到一系列时间帧 `frames`，每一帧是一个波函数 $\psi_{3,frame}(x)$。

我这里的 survival（更准确说是势阱区域内的概率）定义为：

$$
P_{\text{inwell}}(t)=\int_{x\le x_{top}} |\psi(x,t)|^2\,dx
$$

其中 $x_{top}$ 是倾斜势能的势垒顶位置，我将阱内区域取为 $x\le x_{top}$。这个定义实现直接、数值上较稳定，并且可以与后面基于通量的结果相互对照。

![FIGURE: part_d_survival_vs_time](figures/part_d_survival_vs_time.png)

### 用指数拟合估计逃逸率 $\Gamma_n$

如果在一段时间里逃逸可近似为常数速率过程，则有

$$
P_{\text{inwell}}(t)\approx A\,e^{-\Gamma t}\quad\Rightarrow\quad \ln P_{\text{inwell}}(t)\approx \ln A - \Gamma t
$$

我没有从 $t=0$ 开始拟合，而是取了后半段的帧（`last_n = int(len(frames)*0.5)`）。原因是初态并非倾斜势的本征态，演化初期通常存在一个非指数的过渡段。对后半段数据，我用 `scipy.stats.linregress` 对 $\ln P_{\text{inwell}}$ 与 $t$ 做线性回归，斜率的相反数作为逃逸率 $\Gamma$ 的估计。

本次拟合的结果显示：$\Gamma \approx 0.0002201422482195$。

![FIGURE: part_d_log_survival_linear_fit](figures/part_d_log_survival_linear_fit.png)

### 概率通量：势垒顶处的流出强度

除了积分 survival，我还想直接看流出：对一维波函数，概率流密度是

$$
j(x,t)=\frac{\hbar}{m}\,\mathrm{Im}\left(\psi^*(x,t)\,\partial_x\psi(x,t)\right)
$$

我在实现里用 `np.gradient` 计算空间导数，得到 $j(x)$ 的离散数组（见 [`probability_flux`](quantum_tunneling/tdse.py) / [`probability_current`](quantum_tunneling/observables.py)）。随后取势垒顶所在的网格点 $x_{top}$，得到随时间变化的通量序列 $j(x_{top},t)$。

为了与 survival 建立对应关系，我画了 $j(x_{top},t) / P_{\text{inwell}}(t)$。在指数衰减近似成立的稳定区间内，该比值应接近常数，可作为对 $\Gamma$ 的另一种交叉检验。



![FIGURE: part_d_flux_over_survival](figures/part_d_flux_over_survival.png)

### 波函数时间演化动画：直观检查逃逸过程

仅凭一条 survival 曲线，很难区分真实的外泄与数值误差（例如范数损失、边界反射等）。因此我制作了时间演化动画，将倾斜势能、能级线，以及 $\psi$ 的实部、虚部与概率密度叠加在同一张图上，以便直观检查演化过程。

这段动画用 `matplotlib.animation.FuncAnimation` 输出为 gif（notebook 里 `ani.save(...)` 那一段）。主要关注两点：

- 波函数在势垒外侧是否出现持续的外泄，并在 CAP 区域被有效吸收
- CAP 是否主要作用于边界区域，且对阱内波形不产生明显的非物理扰动

![FIGURE: part_d_time_evolution_animation](figures/part_d_time_evolution_animation.gif)

### 偏离 WKB：画 $\ln\Gamma_{TDSE}$ vs $S(E)$

最后我将数值得到的逃逸率与 WKB 作用量联系起来。在最简的 WKB 近似下，有

$$
\Gamma \propto e^{-2S/\hbar}\quad\Rightarrow\quad \ln\Gamma \approx -\frac{2}{\hbar}S + \text{const}
$$

我的做法是：对所有仍位于势垒顶以下的本征态（`i in range(highest_bound_state)`），分别运行 TDSE 得到 $\Gamma_i$；然后在同一外场 $F$ 下，计算对应的 $S(E_i)$。最后将 $(S(E_i),\ln\Gamma_i)$ 作图，并对 $\ln\Gamma$ 与 $S$ 做线性拟合，用来检验是否近似满足线性关系。

![FIGURE: part_d_lnGamma_vs_S](figures/part_d_lnGamma_vs_S.png)
 
本次拟合的结果显示：$ln\Gamma_n=-0.0903*S(E_n)-7.9411$,r 为 -0.86, 有较强的负相关性, 斜率约为 -0.0903，与理论预期 $-2/\hbar = -2$ 有较大偏差。