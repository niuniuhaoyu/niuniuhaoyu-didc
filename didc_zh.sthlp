{smcl}
{* *! version 0.2.1  30sep2026}{...}
{hline}
{p 4 8 2}{bf:didc} —— 差中差断点（Difference-in-Discontinuities）{right:版本 0.2.1}
{hline}

{p 4 4 2}{it:英文帮助：} {help didc}、{help didc_test}、{help didc_bounds}

{title:标题}

{p 4 4 2}
{cmd:didc} —— 当关注的处理与一个混淆政策由同一阈值赋值时，估计处理效应（DiDC 设计），
并作推断；配套 {help didc_zh} 内的 {cmd:didc_test}（效度检验）与 {cmd:didc_bounds}
（部分识别）。方法出自 Picchetti、Pinto & Shinoki (2026, arXiv:2405.18531)。

{title:didc 语法}

{p 8 12 2}
{cmd:didc} {it:depvar} {ifin} {cmd:,} {cmd:runvar(}{it:varname}{cmd:)}
{cmd:time(}{it:varname}{cmd:)} {cmd:pre(}{it:numlist}{cmd:)} {cmd:post(}{it:numlist}{cmd:)}
{cmdab:design(}{it:panel}{cmd:|}{it:rcs}{cmd:)} {cmdab:id(}{it:varname}{cmd:)}
{cmdab:cutoff(}{it:#}{cmd:)} {cmdab:p(}{it:#}{cmd:)} {cmdab:q(}{it:#}{cmd:)}
{cmdab:kernel(}{it:kern}{cmd:)} {cmdab:bwselect(}{it:rule}{cmd:)} {cmdab:h(}{it:#}{cmd:)}
{cmdab:b(}{it:#}{cmd:)} {cmdab:level(}{it:#}{cmd:)} {cmdab:graph(}{it:name}{cmd:)}
{cmd:nolemma1} {cmdab:engine(}{it:rdrobust}{cmd:|}{it:mata}{cmd:)}}

{title:didc 描述}

{pstd}
估计量是{hi:结果随时间之差}在断点两侧的局部多项式回归，即论文引理 1 的左边；右边是两期
断点之差。二者在总体中相等，但有限样本下只有使用共同带宽才相等——{cmd:didc} 始终用共同
带宽，使其一致。

{title:didc 主要选项}

{p 4 8 2}{cmd:runvar(}{it:varname}{cmd:)} 运行变量。必须单位内{hi:时间不变}；命令会检查
并报告违约单位。
{p 4 8 2}{cmd:time(}{it:varname}{cmd:)}、{cmd:pre(}{it:numlist}{cmd:)}、
{cmd:post(}{it:numlist}{cmd:)} 时期标识；数据只允许含这两期。
{p 4 8 2}{cmd:design(panel)}（默认）要求同一批单位在两期都出现，需 {cmd:id()}；单位内构造
结果变化，对其跑一次 RD。{cmd:design(rcs)} 用于重复截面，不允许 {cmd:id()}，对每期各跑一次
RD（共同带宽）再作差。
{p 4 8 2}{cmd:cutoff(}{it:#}{cmd:)} 处理状态变化的运行变量取值，默认 0（内部中心化）。
{p 4 8 2}{cmd:p(#)}、{cmd:q(#)} 估计与偏误矫正的局部多项式阶数（{cmd:q()} 须大于
{cmd:p()}），默认 1、2。
{p 4 8 2}{cmd:kernel()}、{cmd:bwselect()}、{cmd:h()}、{cmd:b()} 透传给
{help rdrobust}。
{p 4 8 2}{cmd:graph(}{it:name}{cmd:)} 画出差分结果的 RD 图并存为 {it:name}（仅
{cmd:design(panel)}）。
{p 4 8 2}{cmd:nolemma1} 跳过引理 1 自检。
{p 4 8 2}{cmd:engine(rdrobust|mata)} 选择引擎。默认 {cmd:engine(rdrobust)} 把局部多项式
拟合、偏误矫正与方差委托给 {help rdrobust}；{cmd:engine(mata)} 使用 Mata 自包含内核，点在
{cmd:p(1)}/{cmd:p(2)} 下与 {cmd:rdrobust} 一致到机器精度，方差为精确有限样本方差（稳健标准
误约为 {cmd:rdrobust} 的 0.94 倍）。默认仍用 {cmd:rdrobust}，因为它与已发表输出逐位一致。

{title:didc 存储结果}

{p 4 4 2}
除通常的 {cmd:rdrobust} 结果外：{cmd:e(tau_didc)}（常规估计）、{cmd:e(tau_didc_bc)}
（偏误矫正）、{cmd:e(se_didc)} 与 {cmd:e(ci_didc_l)}/{cmd:e(ci_didc_r)}（稳健 SE 与区间）、
{cmd:e(mu_plus)}/{cmd:e(mu_minus)}（断点两侧截距）、{cmd:e(lemma1_ok)}（两参数化是否一致）、
{cmd:e(engine)}、以及带宽 {cmd:e(bw_h_l)}、{cmd:e(bw_h_r)}、{cmd:e(bw_b_l)}、
{cmd:e(bw_b_r)}（{warning:不要找 {cmd:e(b)}，它存的是估计值}）。

{title:didc_test 语法}

{p 8 12 2}
{cmd:didc_test} {it:depvar} {ifin} {cmd:,} {cmd:runvar(}{it:varname}{cmd:)}
{cmd:time(}{it:varname}{cmd:)} {cmdab:test(}{it:wald}{cmd:|}{it:ks}{cmd:)}
{cmd:pre(}{it:numlist}{cmd:)} {cmdab:cutoff(}{it:#}{cmd:)} {cmdab:p(}{it:#}{cmd:)}
{cmdab:kernel(}{it:kern}{cmd:)} {cmdab:bw(}{it:#}{cmd:)} {cmdab:degree(}{it:#}{cmd:)}
{cmdab:reps(}{it:#}{cmd:)} {cmdab:seed(}{it:#}{cmd:)} {cmdab:level(}{it:#}{cmd:)}
{cmd:noisily}

{title:didc_test 描述}

{pstd}
{cmd:test(wald)} 检验断点处的跳跃在每个处理前期是否相同（堆叠局部线性 RD 上的联合 Wald
检验）。{cmd:test(ks)} 检验最近两个处理前期的条件均值 E[Y|Z] 是否相同（多项式基上的两样本
KS，bootstrap 施加零假设）。任一被拒绝都警示假设 4（混淆效应时不变）不成立——这是整个设计
所依赖的假设。{cmd:pre()} 须含{bf:至少两期}处理前期。

{title:didc_test 存储结果}

{p 4 4 2}
{cmd:test(wald)}：{cmd:r(F)}/{cmd:r(df)}/{cmd:r(p)}（水平跳跃相等）、
{cmd:r(F_theta)}/{cmd:r(p_theta)}（斜率变化）、{cmd:r(F_joint)}/{cmd:r(p_joint)}（联合）、
{cmd:r(gamma)}/{cmd:r(theta)}、{cmd:r(bw)}/{cmd:r(N)}。
{p 4 4 2}
{cmd:test(ks)}：{cmd:r(ks_above)}/{cmd:r(ks_below)}、
{cmd:r(p_above)}/{cmd:r(p_below)}、{cmd:r(N_above)}/{cmd:r(N_below)}、
{cmd:r(period_new)}/{cmd:r(period_old)}。

{title:didc_bounds 语法}

{p 8 12 2}
{cmd:didc_bounds} {it:depvar} {ifin} {cmd:,} {cmd:runvar(}{it:varname}{cmd:)}
{cmd:time(}{it:varname}{cmd:)} {cmd:pre(}{it:numlist}{cmd:)} {cmd:post(}{it:numlist}{cmd:)}
{cmdab:design(}{it:panel}{cmd:|}{it:rcs}{cmd:)} {cmdab:id(}{it:varname}{cmd:)}
{cmdab:cutoff(}{it:#}{cmd:)} {cmdab:c1(}{it:numlist}{cmd:)} {cmdab:c2(}{it:numlist}{cmd:)}
{cmdab:p(}{it:#}{cmd:)} {cmdab:q(}{it:#}{cmd:)} {cmdab:kernel(}{it:kern}{cmd:)}
{cmdab:bwselect(}{it:rule}{cmd:)} {cmdab:h(}{it:#}{cmd:)} {cmdab:b(}{it:#}{cmd:)}
{cmdab:level(}{it:#}{cmd:)} {cmdab:matrix(}{it:name}{cmd:)}}

{title:didc_bounds 描述}

{pstd}
把假设 4（混淆效应时不变）替换为有界变动假设（论文假设 6），报告识别集：
{it:tau}{sub:c} ∈ [ max{ {it:mu}{sup:+} − c1 , {it:tau}{sub:didc} − c2 } ,
min{ {it:mu}{sup:+} + c1 , {it:tau}{sub:didc} + c2 } ]。
{bf:c1} 界定被混淆结果的时间趋势；{bf:c2} 界定混淆效应本身的漂移（c2 = 0 即恢复假设 4）。
命令报告整个 {cmd:c1()}×{cmd:c2()} 网格与{hi:破坏值}（breakdown）——在另一参数取零时，
某个参数大到多少就无法对效应定号。

{pstd}
{bf:下界大于上界 = 空识别集，这是结果而非报错。} 在 {cmd:c1}={cmd:c2}=0 时，只要
{it:mu}{sup:-} 不恰好为 0 集合就为空——这等价于 {cmd:didc_test} 拒绝假设 4。实践中原点格
几乎总是空的，有信息量的是破坏值。

{title:didc_bounds 存储结果}

{p 4 4 2}
{cmd:e(bounds_lb)}、{cmd:e(bounds_ub)}（下/上界矩阵，行 = c1，列 = c2）、
{cmd:e(bounds_empty)}、{cmd:e(mu_plus)}/{cmd:e(mu_minus)}/{cmd:e(tau_didc)}、
{cmd:e(breakdown_c1)}、{cmd:e(breakdown_c2)}、{cmd:e(identified_set_empty_at_zero)}。

{title:示例}

{p 8 8 2}{cmd:. didc y, runvar(z) time(t) pre(0) post(1) id(id)}
{p 8 8 2}{cmd:. didc_test y, runvar(z) time(t) test(wald) pre(-1 0)}
{p 8 8 2}{cmd:. didc_test y, runvar(z) time(t) test(ks) pre(-1 0)}
{p 8 8 2}{cmd:. didc_bounds y, runvar(z) time(t) pre(0) post(1) id(id) c1(0(0.5)3) c2(0(0.5)3)}

{title:参考文献}

{pstd}
Calonico, S., M. D. Cattaneo, and R. Titiunik (2014). {it:Econometrica} 82(6): 2295-2326.
Grembi, V., T. Nannicini, and U. Troiano (2016). {it:AEJ: Applied Economics} 8(3): 1-30.
Picchetti, P., C. C. X. Pinto, and S. T. Shinoki (2026). {it:arXiv preprint} arXiv:2405.18531.

{title:另见}

{p 4 4 2}
英文帮助：{help didc}、{help didc_test}、{help didc_bounds}、{help rdrobust}
