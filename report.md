# Project 2 研究报告

## 选择

我选了Cusp Potential作为研究对象，因为它在0点有一个尖锐的拐角，这种奇异结构会对波函数的行为产生显著影响。Cusp Potential的数学形式为：
$$
V(x) = V_0 |x|^{\alpha}, \quad 0 < \alpha < 1
$$

## Part A Bound States and Wavefunction Structure

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

结果比较有趣：能级越高，禁阻区概率反而越小。我想这是可能因为Vcusp在远离中心时增长得很快，导致高能态的波函数主要集中在势阱内，难以渗透出去。

![FIGURE/part_a_forbidden_region_probability_vs_n](figures/part_a_forbidden_region_probability_vs_n.png)
![FIGURE/part_a_forbidden_vs_n](figures/part_a_forbidden_vs_n.png)

### 波函数伸出来多远

定义一个"穿透深度"：从转折点开始，波函数的概率密度衰减到 1/e 时的距离。计算这个深度对不同 n 的变化。

高阶态的穿透深度更大，也就是说波函数在禁阻区衰减得更缓。对应的，衰减速率（粗略定义为幅度下降的快慢）也就更小。

![FIGURE: part_a_penetration_decay_vs_n]()




