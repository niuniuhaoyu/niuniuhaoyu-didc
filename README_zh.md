# didc

[English](README.md) | [简体中文](README_zh.md)

**面向 Stata 的差中差断点（Difference-in-Discontinuities）**

[![Stata 17+](https://img.shields.io/badge/Stata-17%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 0.1.0](https://img.shields.io/badge/Version-0.1.0-green.svg)](CHANGELOG.md)
[![Requires: rdrobust](https://img.shields.io/badge/requires-rdrobust-informational.svg)](https://rdpackages.github.io/rdrobust/)

`didc` 为 Stata 提供**差中差断点（DiDC）**的估计、推断、效度检验与部分识别，遵循
Picchetti、Pinto & Shinoki (2026)《Difference-in-Discontinuities: Estimation,
Inference and Validity Tests》（[arXiv:2405.18531](https://arxiv.org/abs/2405.18531)）。

---

## 它解决的问题

在许多政策情境中，关注的处理由运行变量上的一个阈值决定——**同一个阈值还决定了另一件事**。
意大利 5,000 人口以下的市镇在 2001 年放松了财政规则，同时也有更高的市长工资，二者都由
同一个"5,000 人口"阈值决定。其他例子：同一条地理边界上叠加了多个项目；学校阈值同时触发
班级规模规则与经费规则。

此时标准断点回归（RD）失效：运行变量在断点处的**均值潜在结果**不连续，因为*混淆变量*
也在那里跳变。双重差分通常也用不了，因为阈值两侧的单位在水平上不可比。

DiDC 用**时间上的变化**解决它：一个随时间不变的混淆效应，在处理前与处理后贡献同样的
跳跃，在差分中被消去。

## 安装

```stata
* 依赖
ssc install rdrobust, replace
```

然后安装 `didc` 本身。想要结果可复现就用**固定 release**，想用最新就用开发版：

```stata
* 固定 release（推荐）
net install didc, from("https://raw.githubusercontent.com/niuniuhaoyu/niuniuhaoyu-didc/v0.2.1/") replace

* 开发版
net install didc, from("https://raw.githubusercontent.com/niuniuhaoyu/niuniuhaoyu-didc/main/") replace
```

`didc` 需要明确拒绝一种常见失败：若 `rdrobust` 缺失，它会明说并告诉你如何安装，
而不是报 `r(111)`。

用 `which didc`、`which didc_test`、`which didc_bounds` 检查安装，然后运行 `help didc`。

### 内置引擎

`didc` 提供两个估计引擎。`engine(rdrobust)`（默认）把局部多项式估计与推断委托给
`rdrobust`，与原文一致。`engine(mata)` 使用 Mata 中自包含的内核：

```stata
didc y, runvar(z) time(t) pre(0) post(1) id(id) engine(mata) h(0.4) b(0.6)
```

给定 `h()` 与 `b()` 时，`engine(mata)` 完全不调用 `rdrobust`；不给定则仅用 `rdrobust`
做**带宽选择**。

| | `engine(rdrobust)` | `engine(mata)` |
|---|---|---|
| 点估计与偏误矫正 | 参考实现 | `p(1)` 与 `p(2)` 下与之一致到 **1e-16**，已对三角、均匀、Epanechnikov 核验证 |
| 方差 | CCT 的渐近方差，`sigma^2` 由近邻估计 | 同一线性泛函的**精确有限样本方差** |
| 稳健标准误，模型 1 | 参考实现 | 为 `rdrobust` 的 0.94 倍 |
| 覆盖率，模型 1，`n=1000`，150 次重复 | 0.960 | **0.953** |
| 带宽选择 | `rdrobust` | `rdrobust`，或用户给定的 `h()`、`b()` |
| 多项式阶数 | 任意 `p` | 任意 `p` |

两个引擎在关键处数值等价：内置引擎把每个点估计复现到机器精度——包括 `p(2)` 下，两引擎
在六个量上同时一致到 `8.3e-16`——其区间覆盖率 0.953 对 `rdrobust` 的 0.960，差异
完全在 Monte Carlo 误差内。标准误本身残留的 0.94 差异，来自两套 `sigma^2` 约定——CCT
用近邻估计，内置引擎用拟合自身的残差。`engine(rdrobust)` 仍是默认，因为它是参考实现，
且与已发表的 `rdrobust` 输出逐位一致——这正是核查你数字的读者会去跑的东西。

完全弥合标准误差距、以及写一个带宽选择器，是接下来的工作；见
`docs/specs/2026-09-30-didc-v2-design.md`。

## 快速上手

```stata
use "didc_sim1.dta", clear          // 包内自带的模拟数据
didc y, runvar(z) time(t) pre(0) post(1) id(id)

didc_test y, runvar(z) time(t) pre(-1 0) test(wald)   // 混淆效应是否稳定？
didc_test y, runvar(z) time(t) pre(-1 0) test(ks)     // E[Y|Z] 是否稳定？

didc_bounds y, runvar(z) time(t) pre(0) post(1) id(id) c1(0(0.5)5) c2(0(0.5)5)
```

完整、可复现的演练见 `examples/didc_example.do`。

## 数据布局

命令接受**长格式**数据：每单位每期一行。

| `design()` | 含义 | `id()` | 估计什么 |
|---|---|---|---|
| `panel`（默认） | 同一批单位出现在两期 | 必填 | 单位内构造 `DeltaY = y(post) - y(pre)`，再对 `(DeltaY, Z)` 跑**一次** RD |
| `rcs` | 重复截面；两期是独立样本 | 不允许 | 对每期各跑**一次** RD（共同带宽），再作差 |

`runvar()` 必须单位内**时间不变**；命令会检查并告诉你哪些单位违约。

## 方法

### 目标参数

令 `Z` 为运行变量，`z0 = cutoff()`，`D_{i,0}` 为混淆政策，`D_{i,1}` 为关注的处理，二者
都由 `1{Z >= z0}` 赋值。目标参数为

```
tau_c = E[ Y_{i,1}(1,1) - Y_{i,1}(1,0) | Z_i = z0 ]
```

即在断点处、对暴露于混淆变量的单位而言，关注处理的效应。在论文的假设 1–4 下，它等于
DiDC 目标参数

```
tau_didc = tau_1^RD - tau_0^RD
```

### 引理 1：估计量等于对差分做一次 RD

论文的关键计算结果是：*两个 RD 之差*等于*差分的 RD*：

```
tau_didc = lim_{z->z0+} E[ dY | Z=z ] - lim_{z->z0-} E[ dY | Z=z ],   dY = y_1 - y_0
```

所以估计量就是对 `dY` 关于 `Z` 的局部多项式拟合，与 Cattaneo、Calonico & Titiunik (2014)
完全一样——也就是 `rdrobust` 所做的事。

**一个 `didc` 替你处理的有限样本细节。** 该恒等式是关于总体极限的陈述。如果你分别估计
两边、让各自回归选择自己的 MSE 最优带宽，两种参数化在数值上*并不*一致。在模拟数据上：

| | 带宽 | 估计 |
|---|---|---|
| `dY` 的 RD | 0.2165 | 0.300031 |
| `RD(y_post) - RD(y_pre)` | 各自选 | 0.290737 |
| 两者，共同带宽 | 0.2165 | 0.300031 vs 0.300031 |

因此 `didc` 总是**用共同带宽**报告两种参数化，并返回 `e(lemma1_ok)` 作为内置检查。
没有其他 Stata 包这样做。

### 推断

局部多项式估计，配 MSE 最优插件带宽选择器、偏误矫正、以及 Calonico、Cattaneo & Titiunik
(2014) 的稳健置信区间。这些都来自 `rdrobust`；见下方*实现说明*。

### 效度检验（`didc_test`）

| 检验 | 问题 | 设定 |
|---|---|---|
| `test(wald)` | 每个处理前期的断点是否相同？ | 对 `pre()` 中的各期做堆叠局部线性 RD，含各期水平跳跃 `gamma_k` 与斜率变化 `theta_k`；联合 Wald 检验 |
| `test(ks)` | 最近两个处理前期的 `E[Y\|Z]` 是否相同？ | 多项式基上的两样本 Kolmogorov-Smirnov 统计量，bootstrap **施加零假设** |

任一检验被拒绝，都会动摇假设 4（混淆效应的时不变性），从而动摇 DiDC 设计本身。
这两个检验在 Stata 生态里都不存在。

**`test(wald)` 内部一个重要的方差修正。** 论文用核加权的加权最小二乘估计堆叠回归，并用
通常的 `F` 统计量检验。把核权重当 `aweights` 交给 Stata，会让 Stata 把它们当作逆方差
权重，从而低估方差约 `mean(w)` 倍，使检验过于宽松——在论文的模型 1 上实测 5% 名义水平
下拒绝率约 **14%**。因此 `didc_test` 用 `matrix accum` 构造正确的同方差 WLS 方差
`sigma2 (X'WX)^{-1}`，其中 `sigma2 = sum(w u^2)/(sum(w) - k)`，并直接构造 Wald 统计量。
这恢复了名义水平，是朝论文*所述*检验方向、偏离论文*实现*的一处改动。

### 部分识别（`didc_bounds`）

当假设 4 不可信时，用有界变动（假设 6）替代它，报告识别集而不是一个点：

```
tau_c in [ max{ mu_plus - c1 , tau_didc - c2 } ,
           min{ mu_plus + c1 , tau_didc + c2 } ]
```

其中 `mu_plus = dY^+`。`c1` 界定被混淆结果中的时间趋势；`c2` 界定混淆效应本身的漂移。
命令报告整个 `c1 x c2` 网格与**破坏值（breakdown）**：`c1` 或 `c2` 可以大到多少，
仍能对效应定号。

**空识别集是特性，不是 bug。** 在 `c1 = c2 = 0` 时，只要 `dY^- != 0`（真实数据中几乎总是
如此），集合就是空的。这是反对假设 4 的直接证据，`didc_bounds` 会这样告诉你，而不是报错。

## 实现说明（哪些是原创、哪些不是）

| 组件 | 如何产生 | 原创？ |
|---|---|---|
| panel/rcs 处理、`dY` 构造、时不变性检查 | 本项目编写 | 是 |
| 点估计、MSE 最优带宽、偏误矫正、稳健 CI | `engine(rdrobust)` 委托 `rdrobust`；`engine(mata)` 自行计算 | 否 / **是** |
| 引理 1 自检 | 本项目编写 | 是 |
| 堆叠 RD 的 Wald 检验 | 本项目编写 | 是，Stata 无先例 |
| bootstrap KS 检验 | 本项目编写 | 是，Stata 无先例 |
| 部分识别、破坏值、空集检测 | 本项目编写 | 是，Stata 无先例 |

**内置 Mata 引擎**（`_didc_mata.ado`，自 v0.2.0 起）在不调用 `rdrobust` 的情况下计算局部
多项式拟合、偏误矫正与方差。其点估计对 `p(1)`、三种核都复现 `rdrobust` 到 `1e-16`，
即 `examples/_test_engine.do` 中的检查。带宽选择仍委托，`p(2)` 及以上被拒绝而非近似；
见上表。

委托估计内核正是原文的做法——它明确建立在 Calonico、Cattaneo & Titiunik (2014) 之上。
手搓带宽选择器与稳健方差会产生无法与已发表数字核对的差异，因此 mata 引擎与 `rdrobust`
路径并存，并报告实测差距。

自包含估计（带宽选择器、局部多项式、偏误矫正与稳健方差全部用 Mata 实现）是 v2 的路线；
落地后两套实现互相校验。

### 与 `RDDID` 的区别

`RDDID`（SSC `s459586`，Jonathan Dries）做的是*横截面* DiDC：处理组的断点减去对照组的
断点。它没有时间维度、没有效度检验、没有部分识别。`didc` 估计的是 Grembi 等 (2016) 的
*时间差分*设计，是另一种研究设计。两者只共享 `rdrobust` 作为计算后端，再无其他。

### 原文中的不一致及 `didc` 的处理

完整记录见 `docs/research-notes.md` 第 9 节。实质性几处：

1. **公式 (7)** 漏印了水平项 `gamma_k D`，那样只留下断点处的斜率变化、没有对断点本身的
   检验。`didc_test` 同时估计两者并都报告；头条检验针对水平跳跃。
2. **`KS_-`** 印成对 `Z_i >= z0` 求和，与其自身"断点以下支撑"的定义矛盾。`didc_test`
   改用 `Z_i < z0`，与 `KS_+` 对称。
3. **引理 6** 被标注在"假设 1、2、5"下；其推导需要假设 3 与假设 6。`didc_bounds` 遵循
   推导。
4. **表 1–4** 标题写"10,000 次 Monte Carlo 实验"，而第 6 节说 1,000。
   `examples/_test_mc_coverage.do` 把重复次数作为参数，两种都能跑。

## 存储结果

`didc` 返回通常的 `rdrobust` 结果，外加：

| `e()` | 内容 |
|---|---|
| `e(tau_didc)` | 常规 DiDC 估计 |
| `e(tau_didc_bc)` | 偏误矫正估计 |
| `e(se_didc)`、`e(ci_didc_l)`、`e(ci_didc_r)` | 稳健标准误与区间 |
| `e(mu_plus)`、`e(mu_minus)` | `dY^+` 与 `dY^-`，左/右截距 |
| `e(lemma1_ok)` | 若两种参数化在共同带宽下一致则为 1 |
| `e(tau_delta_of_rds)` | "RD 之差"参数化 |
| `e(design)`、`e(cutoff)`、`e(n_units)`、`e(n_below)`、`e(n_above)` | 设计与样本 |
| `e(bw_h_l)`、`e(bw_h_r)`、`e(bw_b_l)`、`e(bw_b_r)` | 带宽（注意：**不是** `e(b)`，后者存估计值） |

`didc_test` 对 `test(wald)` 返回 `r(p)`、`r(F)`、`r(df)`、`gamma` 与 `theta` 矩阵以及
斜率/联合检验；对 `test(ks)` 返回 `r(ks_above)`、`r(ks_below)`、`r(p_above)`、`r(p_below)`。

`didc_bounds` 返回 `e(bounds_lb)`、`e(bounds_ub)`、`e(bounds_empty)`、
`e(breakdown_c1)`、`e(breakdown_c2)`。

## 数字正确的证据

| 检查 | 文件 | 结果 |
|---|---|---|
| 线性差分结果上的精确恢复与带宽不变性 | `examples/_test_closedform.do` | 精确到 `1e-10` |
| 混淆变量被消去：`didc` 恢复 `tau`，而仅用处理后的 RD 返回 `c + tau` | `examples/_test_confounder.do` | 与预期一致 |
| **独立**局部线性估计（经 `regress` 的加权最小二乘，不用 `rdrobust`） | `examples/_test_vs_wls.do` | 两种核都一致到 `0.0000000000` |
| 经 `reshape`（而非包自身路径）重建 `dY` | `examples/_test_vs_wls.do` | 完全一致 |
| Monte Carlo 偏误与覆盖率，模型 1–4 | `examples/_test_mc_coverage.do` | 见下 |
| 两个效度检验的 size 与 power | `examples/_test_validity_size.do` | KS：size `0.050`，power `1.000`；Wald：size `0.025`，power `0.990` |
| 闭式界、空集、破坏值 | `examples/_test_bounds.do` | 逐格精确 |

Monte Carlo，200 次重复，`n = 1000`，`tau = 0`，对照论文表 1–4。"RD"指对处理后横截面
直接跑的普通断点回归，即忽略混淆变量：

| 模型 | didc 偏误 | RMSE | 覆盖率 | RD 偏误 | RD 覆盖率 | 论文：didc 偏误 / 覆盖率 |
|---|---|---|---|---|---|---|
| 1，有混淆，固定形式 | -0.004 | 0.052 | **0.945** | 0.541 | **0.000** | -0.002 / 0.944 |
| 2，无混淆，固定形式 | -0.004 | 0.052 | 0.945 | 0.041 | 0.860 | -0.003 / 0.944 |
| 3，有混淆，时变形式 | -0.008 | 0.052 | **0.950** | 0.541 | **0.000** | -0.001 / 0.937 |
| 4，无混淆，时变形式 | -0.008 | 0.052 | 0.950 | 0.041 | 0.860 | -0.008 / 0.937 |

论文表 2 对模型 2 报告 RD 偏误 0.043、覆盖率 0.861；同样设计此处给出 0.041 与 0.860。
模型 1 和 3 的 RD 偏误此处更小（0.541 对 1.043），因为论文未报告其模拟中混淆变量的
量级，而本包取 0.5。定性结论被精确复现：当断点处存在混淆变量时，普通 RD 严重偏误、
区间永不稳定覆盖；而 `didc` 居中，且即使函数形式随时间变化，覆盖率仍接近名义水平。

一次运行整套：

```stata
do "examples/_test_all.do"        // 或： do "examples/_test_all.do" 1000
```

`examples/reference/run_rdrobust.R` 在有 R 的机器上，用 **R** 包 `rdrobust` 复现
`_test_vs_wls.do` 的对比；仓库内的检查不需要 R。

## 引用

方法的 BibTeX：

```bibtex
@article{picchetti2026didc,
  author  = {Picchetti, Pedro and Pinto, Cristine C. X. and Shinoki, St{\'e}phanie T.},
  title   = {Difference-in-Discontinuities: Estimation, Inference and Validity Tests},
  journal = {arXiv preprint arXiv:2405.18531},
  year    = {2026}
}
```

实证例子遵循：

```bibtex
@article{grembi2016do,
  author  = {Grembi, Veronica and Nannicini, Tommaso and Troiano, Ugo},
  title   = {Do Fiscal Rules Matter?},
  journal = {American Economic Journal: Applied Economics},
  volume  = {8},
  number  = {3},
  pages   = {1--30},
  year    = {2016}
}
```

## 许可

AGPL-3.0。见 `LICENSE`。
