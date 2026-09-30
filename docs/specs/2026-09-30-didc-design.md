# didc — Stata 包设计（Spec）

> 日期：2026-09-30
> 状态：待用户审阅（Draft）
> 作者：Haoyu Niu
> 定位：Difference-in-Discontinuities（DiDC，断点处存在混淆政策时的双重差分断点设计）的 Stata 实现
> 前置研究：`didc/docs/research-notes.md`（估计量精确定义、假设、DGP、裁决留痕）

---

## 1. 背景与目标

**方法**：Picchetti, Pinto & Shinoki, *Difference-in-Discontinuities: Estimation, Inference and
Validity Tests*，arXiv:2405.18531（v2, 2026-01）。

**问题**：当阈值处除了关注政策之外还**同时存在一个混淆政策**（例如意大利市政：市长薪资与财政规则放松共用 5000 人口这个门槛），标准 RDD 的连续性假设失效，估计被系统污染；而 DiD 又不适用（阈值两侧单位本身不可比、平行趋势不成立）。DiDC 用时差消掉**时不变的混淆跳跃**，恢复无偏估计。

**缺口（已核实）**：
- SSC 全量目录 3402 条逐条比对：**无 DiDC 主题命令**；`difference-in-discontinuities` 仅 3 处命中，全部指向 `RDDID`（S459586）。
- Stata `RDDID`（Jonathan Dries）是**横截面「处理组断点 − 对照组断点」版**的薄封装（描述原文：*"by wrapping the rdrobust package"*），**不含**：时间差分设计、CCT 偏误校正与稳健 CI、混淆效应时不变检验、KS 函数形式检验、部分识别与敏感度分析。
- 论文本身未发布 R/Stata 包（其估计直接调用 `rdrobust`）。

**目标**：实现一个**正确、可用、有文档、可复现**的 Stata 包 `didc`，把上述空洞一次补齐。

**差异化对齐（本包的存在理由）**：论文把估计量委托给 `rdrobust` 这一事实，给了本包**最硬的验证锚点**——`didc` 的点估计与稳健 CI 必须与「对 `(ΔY, Z)` 直接调用 R `rdrobust`」逐位对齐（容差 ≤ 1e-6）。这与 `contdid` 对照 R `contdid` 的策略同源，是本包相对生态中其他「快速移植包」的核心竞争力。

---

## 2. 方法概述（高层）

1. 计算每个单位的结果变化 `ΔY_i = Y_{i,1} − Y_{i,0}`。
2. 在阈值 `z_0` 两侧分别对 `(ΔY, Z)` 做 **p 阶局部多项式**，取 `z_0` 处左右极限之差：

   `τ^DiDC = Δμ_+ − Δμ_−`

   依据 `Lemma 1`（「RDD 之差」= 「差分之 RDD」）：**只需一次 RDD**，不必分两期估。
3. 带宽用 CCT 的 MSE-最优 plug-in 选择；推断用 **CCT 偏误校正 + 稳健置信区间**（Theorem 1）。
4. 两条**有效性检验**：混淆效应是否时不变（stacked-RD Wald）；潜在结果函数形式是否时不变（两样本 KS，recentered bootstrap）。
5. 假设 A4 不成立时，给出**部分识别边界**（有界变异）与 breakdown value。

---

## 3. v1 范围

### 做（v1）

| 模块 | 内容 |
|---|---|
| 估计 | 两期 DiDC：`ΔY` 上的局部多项式（`ν = 0`，`p` 可选 1–3），CCT 偏误校正与稳健 CI |
| 带宽 | `mserd` / `msetwo` / 手动 `h()`、`b()`；复用 `rdrobust` 的 plug-in 选择器 |
| 检验一 | Stacked-RD：`H_0: θ_0 = θ_{−1} = … = θ_{−K}`，Wald F 检验 |
| 检验二 | 两样本 KS：函数形式时不变性，recentered bootstrap 出的 p 值 |
| 部分识别 | 有界变异（`Lemma 6`）的 `[τ_c^LB, τ_c^UB]` + breakdown value + `c_1 × c_2` 网格 |
| 输出 | `e()`/`r()` 结果、可选图（RD 图 + ΔY 断点图）、可导出的边界表 |
| 配套 | 模拟数据生成器、英文 README、三份 `sthlp`、LICENSE、CHANGELOG、`docs/plans`+`docs/specs` |
| 验证 | 闭式单测 + 对照 R `rdrobust`（≤1e-6）+ 复现论文 Table 1–4 + 检验的 size/power |

### 不做（v2+，明确列为未来）

- **模糊 DiDC**（fuzzy：`D` 在阈值处不完全遵从）
- **模性部分识别**（`Lemma 7`–`Lemma 10`：互补/替代 + 有界结果）→ v2，因为需要 `y^min`/`y^max` 与 `Y_1^±`，接口设计要单独打磨
- **A4′ 乘性效应估计量**（`Lemma 3`；论文自述 future work）
- 协变量调整、聚类推断、多期/交错 DiDC
- 二维边界 DiDC（`rd2d` 覆盖另一设计，不重复）

---

## 4. 命令语法（设计）

三个命令，各管一件事（与 `rdrobust` 家族 `rdrobust`/`rddensity`/`rdlocrand` 的组织方式一致）。

### 4.1 `didc` —— 主估计

```
didc depvar [if] [in] , runvar(varname) time(varname) pre(numlist) post(numlist)
     [ design(panel|rcs) id(varname) cutoff(#) p(#) q(#) kernel(string)
       bwselect(string) h(#) b(#) level(#) graph(string) nolemma1 ]
```

| 项 | 说明 |
|---|---|
| `depvar` | 结局变量（**长表**：每个单位每期一行） |
| `runvar(varname)` | 运行变量 `Z`（要求时不变，命令会校验） |
| `time(varname)` | 期别变量 |
| `pre(numlist)` / `post(numlist)` | 处理前 / 处理后的期别取值（各一期） |
| `design(panel\|rcs)` | 默认 `panel`：同一批单位两期都有观测，**必须给 `id()`**；`rcs`：重复截面（两期是独立样本） |
| `id(varname)` | 单位标识；`design(panel)` 时必填 |
| `cutoff(#)` | 阈值 `z_0`，默认 **0**（内部把 `Z` 中心化到 0） |

#### 4.1.1 两种数据结构（**必须区分**）

论文假设面板（同一批 `i` 两期），估计量建立在 `E[ΔY | Z = z]` 上。`ΔY` 的构造方式取决于数据：

| | `design(panel)`（默认） | `design(rcs)` |
|---|---|---|
| 数据 | 同单位两期成对出现 | 两期为独立样本 |
| `id()` | **必填** | 不可用 |
| 估计 | 先算 `ΔY_i` → 对 `(ΔY, Z)` 跑一次 RD | 两期各跑一次 RD（**强制同一组 `h`/`b`**）→ 相减 |
| 依据 | `Lemma 1` 的左式 | `Lemma 1` 的右式 |
| `Lemma 1` 自检 | 两式都算得出 → **可自检** | 只有一式 → 自检自动关（`e(lemma1_ok) = .`） |

> `rcs` 分支不是额外负担：`Lemma 1` 自检本来就要跑 `y_post` 与 `y_pre` 两个 RD，那正是 `rcs` 的估计式。
> 两种设计返回**同一套 `e()` 结构**，`e(design)` 标明用的是哪一种。

### 4.2 `didc_test` —— 有效性检验

```
didc_test depvar [if] [in] , runvar(varname) time(varname) [cutoff(#)]
     test(wald|ks) pre(numlist) [ p(#) q(#) kernel(string) bw(string) reps(#) seed(#) level(#) ]
```

- `test(wald)`：`pre()` 给出 ≥2 个处理前期的期别，做式 (7) 的 stacked RD 联合 Wald 检验。
- `test(ks)`：`pre()` 给出 ≥2 个处理前期，对左右两侧分别做 KS 检验并给 bootstrap p 值。

### 4.3 `didc_bounds` —— 部分识别与敏感度

```
didc_bounds depvar [if] [in] , runvar(varname) time(varname) pre(numlist) post(numlist)
     [ cutoff(#) c1(numlist) c2(numlist) p(#) kernel(string) bwselect(string) graph ]
```

- `c1()` / `c2()` 给网格（默认各 `0(0.5)5`），输出边界表；额外报告 **breakdown value**（下界跨过 0 的最小 `c_2`）。

### 4.4 返回结果（共用）

`e()`（`didc`）：

| 名称 | 内容 |
|---|---|
| `e(b)` | 点估计（偏误校正后） |
| `e(V)` | 稳健方差 |
| `e(tau_conventional)` / `e(tau_bc)` | 常规 / 偏误校正估计 |
| `e(se_conventional)` / `e(se_robust)` | 对应标准误 |
| `e(ci_conventional)` / `e(ci_robust)` | 置信区间（矩阵） |
| `e(h)` / `e(b)`（带宽） | 主带宽 / 偏误带宽 |
| `e(N_h)` | 有效样本量（左右各一） |
| `e(dy_stats)` | `ΔY` 的描述（均值、左右极限） |

**DiDC 专有输出**（本包相对 `rdrobust` 的新增语义）：

| 名称 | 内容 |
|---|---|
| `e(dy_right)` | `ΔY^+ = e(beta_Y_p_r)[1,1]`（阈值上方 `ΔY` 极限；§6 边界用） |
| `e(dy_left)` | `ΔY^− = e(beta_Y_p_l)[1,1]`（阈值下方 `ΔY` 极限） |
| `e(tau_didc)` | `τ^DiDC = e(dy_right) − e(dy_left) ≡ e(tau_cl)`（恒等，用于自检） |
| `e(lemma1_ok)` | `Lemma 1` 自检：`RD of ΔY` 与 `Δ of RDs`（同一 `h`/`b`）是否逐位一致（`0/1`） |
| `e(tau_delta_of_rds)` | 另一种参数化的值（`RD(y_post) − RD(y_pre)`，强制同带宽） |

> ⚠️ **`e(b)` 名冲突**：偏误带宽不得占用 `e(b)`。裁决——带宽存 `e(bw_h_l)` / `e(bw_h_r)` / `e(bw_b_l)` / `e(bw_b_r)`，点估计存 `e(b)`。写进 sthlp。

---

## 5. 估计量与推断（v1）

1. 校验：`time()` 恰好两期（`pre()`/`post()` 各一个取值）；`runvar()` 时不变（同单位跨期唯一）；`runvar()` 在 `z_0` 两侧都有观测；`depvar` 无缺失（或报出丢弃数）。
2. 生成 `dy = y(post) − y(pre)`（一单位一行）。
3. `Z ← Z − z_0`。
4. 对 `(dy, Z)` 做局部多项式：右侧 `Z ≥ 0`、左侧 `Z < 0`，核 `K`，主带宽 `h`。
5. 点估计 `τ̂ = Δμ̂_+ − Δμ̂_−`；偏误校正 `τ̂^{bc}`（用 `q > p`、带宽 `b`）；稳健方差 `V^{bc}`；CI。
6. 带宽：`bwselect(mserd)` 为默认（CCT 的单带宽选择）；`msetwo` 与手动 `h()`/`b()` 可选。
7. `graph(string)`：`rdplot` 风格的 `ΔY` vs `Z` 断点图（左右分箱 + 拟合线 + 阈值竖线）。

**严格边界**：本包**不自造带宽选择器和方差估计器**，直接复用 CCT 的成熟实现。

### 5.1 实现策略与诚实边界（必须写进 README）

| 组件 | v1 怎么来 | 是否本包原创 |
|---|---|---|
| `ΔY` 构造、面板校验、Z 时不变校验 | 自写 | ✅ |
| 局部多项式点估计 `Δμ^±`、偏误校正 `τ^bc`、稳健 CI | **委托 `rdrobust`**（对 `(ΔY, Z)` 调用） | ❌（与论文做法一致） |
| MSE-最优带宽选择器 | **委托 `rdrobust`**（`mserd`/`msetwo`） | ❌ |
| `Lemma 1` 内置自检（两种参数化一致性） | 自写 | ✅ |
| 混淆效应时不变检验（stacked-RD Wald） | 自写 | ✅ **Stata 生态无先例** |
| 函数形式时不变检验（两样本 KS + recentered bootstrap） | 自写 | ✅ **生态无先例** |
| 部分识别边界 + breakdown value + 空集检测 | 自写 | ✅ **生态无先例** |
| 论文内部不一致的裁决与文档 | 自写 | ✅ |

**为什么这样切**：论文本身就把估计委托给 `rdrobust`（§3 "we rely on the methodologies proposed by Calonico et al. 2014"）。自造带宽/方差只会引入**不可验证**的偏差，反而摧毁本包的护城河（§7 的 1e-6）。因此 v1 把原创性全部压在**设计层、检验层、部分识别层**——这三层恰好是 `RDDID` 与生态里其它包**完全没有**的部分。

**与 `RDDID` 的区别（README 必须正面回答）**：`RDDID` 做的是「处理组断点 − 对照组断点」的**横截面两组**设计，无时间维度、无检验、无部分识别；`didc` 做的是**时间上加差分**的 DiDC 设计，并自带两条有效性检验与敏感度分析。两者共享 `rdrobust` 作为计算后端，但**研究设计不同、统计内容不同**。

---

## 6. 检验与部分识别（v1）

见 `research-notes.md` §5、§6。三条实现约束：

1. **Wald 检验**的带宽不确定性：论文 Appendix D.1 的结论是「任何数据驱动的最优带宽都行」，v1 默认用**各期 CCT 带宽的最小值**，并允许 `bw()` 覆盖。
2. **KS 检验**必须用 **recentered bootstrap**（在原假设下重抽），不能用朴素 bootstrap；`reps()` 默认 999，`seed()` 默认可复现（与 `contdid` 同策略）。
3. **边界推断**：检测是否落在 max/min 的 kink 上；kink 时给出明确警告（标准 bootstrap 失效），并同时报告「若各项不等则标准 bootstrap 有效」的判断。
4. **空识别集**：`LB > UB` 是**真实且必须支持**的输出（等价条件 `ΔY⁻ ≠ 0`，符号无关，见 `research-notes.md` §11.4）。`didc_bounds` 不得报错，须显式标记 `empty` 并解释为「数据与该组敏感度参数不相容，即 `A4` 被数据否决」。注意 `(c1, c2) = (0, 0)` 这一格**几乎总是空的**，因此**必须**同时把 breakdown value 报在显眼位置，避免使用者误读。

---

## 7. 验证方案（本包的护城河，按强度排序）

| # | 检验 | 通过标准 |
|---|---|---|
| V1 | **闭式单测**：`ΔY = a + bZ + τ·1{Z≥0} + ε`（每侧严格线性） | `p=1`、任意带宽下 `τ̂` 与真值完全一致（≤1e-10），且与 `h` 无关 |
| V2 | **混淆消除**：`μ_t` 含时不变跳跃 `c` | `didc` 复原 `τ`；同一数据上「仅用处理后一期跑 RDD」得到 `c + τ`；两者之差 ≈ `c` |
| V3 | **对照 R `rdrobust`**：导出 `(ΔY, Z)` 为 CSV，R 端 `rdrobust(y=dy, x=Z, c=0, ...)` | 点估计、稳健 SE、稳健 CI 全部 **≤ 1e-6** |
| V4 | **复现论文 Table 1–4**：Model 1–4，`n=1000`，1000 次（smoke）与 10000 次（终检） | `τ=0` 下平均偏误 < 0.01、覆盖率 ∈ [0.92, 0.96]；且**同时**复现 RDD 在 Model 1 的崩溃（偏误 > 1、覆盖率 ≈ 0）作为对照 |
| V5 | **检验的 size**：时不变混淆的 DGP，`test(wald)` × 1000 次 | 经验拒绝率 ∈ [0.03, 0.07]（名义 5%） |
| V6 | **检验的 power**：Model 3（时变）DGP | 拒绝率随样本量上升，`n=1000` 时 > 0.5 |
| V7 | **KS 检验 size** | 原假设 DGP 下 bootstrap p 值近似均匀（KS 型检验的 size 检验） |
| V8 | **可移植性冒烟**：`.pkg` 内容清单与实际文件一致、三个命令 `which` 可寻、示例一键可跑 | 全绿 |

> V3 是**发布门槛**：达不到 1e-6 不允许打 tag。

---

## 8. 包结构

```
didc/
├── README.md                      # 英文：Overview / Installation / Quick start / Method / Stored results / Citation / License
├── LICENSE                        # AGPL-3.0（与既有包一致）
├── CHANGELOG.md
├── didc.pkg                       # SSC 包元数据（含 require rdrobust）
├── stata.toc
├── didc.ado                       # 主估计
├── didc_test.ado                  # 有效性检验（wald / ks）
├── didc_bounds.ado                # 部分识别与敏感度
├── _didc_prep.ado                 # 数据校验与 ΔY 构造
├── _didc_lpoly.ado                # 局部多项式包装（委托 CCT）
├── _didc_ks.ado                   # KS 统计量与 recentered bootstrap
├── _didc_display.ado              # 结果表输出
├── didc.sthlp
├── didc_test.sthlp
├── didc_bounds.sthlp
├── examples/
│   ├── didc_simdata.do            # 论文 Model 1–4 的 DGP（可设种子）
│   ├── didc_example.do            # 一键端到端
│   ├── _test_closedform.do        # V1
│   ├── _test_confounder.do        # V2
│   ├── _test_vs_rdrobust.do       # V3（Stata 侧）
│   ├── _test_mc_coverage.do       # V4
│   ├── _test_validity_size.do     # V5/V6/V7
│   └── reference/
│       ├── run_rdrobust.R         # V3 的 R 侧
│       └── *.csv                  # 导出的 (ΔY,Z) 与 R 的参考输出
├── data/
│   └── didc_sim.dta
└── docs/
    ├── research-notes.md
    ├── specs/2026-09-30-didc-design.md
    └── plans/2026-09-30-didc-plan.md
```

> **布局裁决**：采用**扁平结构**（`.ado`/`.sthlp` 置于仓库根目录），与作者已发布并验证过的 `contdid` 保持一致，
> 以保证 `net install didc, from("https://raw.githubusercontent.com/niuniuhaoyu/didc/main/")` 的行为可预期。

---

## 9. 交付物与验收标准

- [ ] 三个命令 + 三份 sthlp，`net install` 可装、`which` 可寻
- [ ] `examples/didc_example.do` 一键复现（含图）
- [ ] **V3 对照 R `rdrobust` 全部 ≤ 1e-6**
- [ ] V1、V2 断言通过；V4 的 Monte Carlo 达到论文量级
- [ ] V5–V7 的 size/power 报告写入 `examples/reference/`
- [ ] 英文 README 完整（含 `research-notes.md` §9 的裁决摘要）
- [ ] AGPL-3.0 + CHANGELOG（记录与论文前提列表不一致处的处理）
- [ ] 推送 GitHub `niuniuhaoyu/didc`；`didc.DAT`/SSC 提交材料备好（不在 v1 验收内）

---

## 10. 边界与未来（YAGNI）

- ❌ 不实现模糊 DiDC、模性边界、A4′ 估计量、多期交错、协变量、聚类
- ✅ v2 候选：`Lemma 7`–`Lemma 10` 的模性边界；fuzzy 设计；协变量调整
- ✅ **v2 候选（用户明确要求，2026-09-30 记录）：自造估计核**——在 Mata 中自行实现
  `rdrobust` 目前承担的三件事：(1) CCT 的 MSE-最优带宽 plug-in 选择器（`mserd`/`msetwo`）、
  (2) 局部多项式点估计与偏误校正、(3) 稳健方差与 CI。
  **动机**：摆脱对 `rdrobust` 的运行时依赖，把「估计核」也变成原创贡献。
  **代价（为什么 v1 不做）**：CCT 的带宽选择器需要估计高阶导数、jackknife 方差与规则性条件，
  工作量大一个数量级；且一旦自造，§7 的「对照 R `rdrobust` ≤ 1e-6」将**从构造性保证变成需要攻克的目标**。
  **建议做法**：v1 发布后再动；自造核写完后，把两份实现（自造 vs `rdrobust`）作为**互验**保留在 `examples/reference/` 中——
  届时那两份实现互相验证，本身就是比现在更强的正确性证据。
- ✅ 与 `contdid` 形成「现代 DiD/RD 工具链」的前两件：`contdid`（连续处理 DiD）+ `didc`（DiDC）

---

## 11. 实施顺序

1. 骨架 + `didc.pkg` + `research-notes.md`（已完成）+ 本 spec
2. `docs/plans/2026-09-30-didc-plan.md`（逐任务、带勾选框）
3. V1 闭式单测（先写测试 → 预期失败）
4. `didc` 主估计（委托 CCT），过 V1 → V2 → V3
5. `didc_test`（wald → ks），过 V5–V7
6. `didc_bounds`，过 V4 的边界表
7. Monte Carlo 终检（V4，10000 次）
8. README / sthlp / CHANGELOG / 示例
9. 推送 GitHub
