{smcl}
{* *! version 0.1.0  30sep2026}{...}
{vieweralsosee "didc" "help didc"}{...}
{vieweralsosee "didc_test" "help didc_test"}{...}
{title:Title}

{phang}
{bf:didc_bounds} --- Partial identification and sensitivity analysis for the
difference-in-discontinuities design

{title:Syntax}

{phang}
{cmd:didc_bounds} {it:depvar} {ifin} {cmd:,} {cmd:runvar(}{it:varname}{cmd:)}
{cmd:time(}{it:varname}{cmd:)} {cmd:pre(}{it:numlist}{cmd:)} {cmd:post(}{it:numlist}{cmd:)}
{hang2}{cmdab:design(}{it:panel}{cmd:|}{it:rcs}{cmd:)} {cmdab:id(}{it:varname}{cmd:)}
{cmdab:cutoff(}{it:#}{cmd:)}
{hang2}{cmdab:c1(}{it:numlist}{cmd:)} {cmdab:c2(}{it:numlist}{cmd:)}
{cmdab:p(}{it:#}{cmd:)} {cmdab:q(}{it:#}{cmd:)} {cmdab:kernel(}{it:kern}{cmd:)}
{cmdab:bwselect(}{it:rule}{cmd:)} {cmdab:h(}{it:#}{cmd:)} {cmdab:b(}{it:#}{cmd:)}
{cmdab:level(}{it:#}{cmd:)} {cmdab:matrix(}{it:name}{cmd:)}}

{title:Description}

{pstd}
{cmd:didc_bounds} replaces Assumption 4 of Picchetti, Pinto & Shinoki (2026) --
that the confounding effect is time-invariant -- with a bounded-variation
assumption, and reports the resulting identified set for the treatment effect on
the confounded units, {it:tau}{sub:c}. This follows Section 5.1 of the paper.

{pstd}
The paper's Assumption 6 introduces two sensitivity parameters:

{pmore}
{bf:c1} bounds {it:| E[Y(1,0) at t=1 - Y(1,0) at t=0 | Z = z] |}: a time trend in
the potential outcome of the confounded units.

{pmore}
{bf:c2} bounds the drift in the {hi:confounding effect} itself between the two
periods. Setting c2 to zero restores Assumption 4.

{pstd}
Under bounded variation (Lemma 6), with {it:mu}{sup:+} = {it:dY}{sup:+} the
intercept of the differenced outcome above the cutoff and {it:tau}{sub:didc} =
{it:mu}{sup:+} - {it:mu}{sup:-},

{pmore}
{it:tau}{sub:c} is in [ max{ {it:mu}{sup:+} - c1 , {it:tau}{sub:didc} - c2 } ,
min{ {it:mu}{sup:+} + c1 , {it:tau}{sub:didc} + c2 } ]

{pstd}
The command reports this set over the whole grid of {cmd:c1()} and {cmd:c2()},
and also reports {hi:breakdown values}: how large each sensitivity parameter can
become, holding the other at zero, before the sign of the effect is no longer
identified.

{title:Options}

{phang}
{cmd:runvar()}, {cmd:time()}, {cmd:pre()}, {cmd:post()}, {cmd:design()},
{cmd:id()}, {cmd:cutoff()}, {cmd:p()}, {cmd:q()}, {cmd:kernel()},
{cmd:bwselect()}, {cmd:h()}, {cmd:b()} and {cmd:level()} are as in
{help didc}. The estimate itself is produced by calling {cmd:didc}.

{phang}
{cmd:c1(}{it:numlist}{cmd:)} is the grid of values for the first sensitivity
parameter; default is {cmd:0(0.5)5}.

{phang}
{cmd:c2(}{it:numlist}{cmd:)} is the grid for the second; default is
{cmd:0(0.5)5}.

{phang}
{cmd:matrix(}{it:name}{cmd:)} stores the lower-bound matrix under {it:name}, with
rows labelled by {cmd:c1()} and columns by {cmd:c2()}.

{title:Reading the output}

{pstd}
{bf:A cell with the lower bound above the upper bound is an EMPTY identified
set, and this is a real result rather than an error.} At {cmd:c1}
= {cmd:c2} = 0 the set collapses to the single point {it:tau}{sub:didc} and is
empty whenever {it:mu}{sup:-} = {it:dY}{sup:-} is not exactly zero: the two
restrictions implied by the assumptions are then mutually inconsistent. That is
the same information as a rejection of Assumption 4 by the Wald test in
{help didc_test}. In practice the origin cell is empty in almost every real
dataset, so the informative quantities are the {hi:breakdown values}, not the
corner cell.

{pstd}
{bf:The bounds use the conventional, not the bias-corrected, estimate.} Lemma 6
plugs the estimator of equation (4) of the paper into the max and min, and that
estimator is the conventional local polynomial one. In Stata terms, the command
uses {cmd:e(tau_didc)} and {cmd:e(mu_plus)}, not {cmd:e(tau_didc_bc)}.

{title:Stored results}

{p2colset 11 20 22 2}{...}
{p2col:{bf:e(bounds_lb)}:matrix of lower bounds; rows are c1, columns are c2}
{p2col:{bf:e(bounds_ub)}:matrix of upper bounds}
{p2col:{bf:e(bounds_empty)}:matrix of 0/1 flags, 1 where the lower bound exceeds the upper bound}
{p2col:{bf:e(mu_plus)}, {bf:e(mu_minus)}, {bf:e(tau_didc)}:the ingredients}
{p2col:{bf:e(breakdown_c1)}:largest c1, at c2 = 0, for which the lower bound stays above zero}
{p2col:{bf:e(breakdown_c2)}:largest c2, at c1 = 0, for which the lower bound stays above zero}
{p2col:{bf:e(identified_set_empty_at_zero)}:1 if the set is empty at c1 = c2 = 0}
{p2col}

{title:Example}

{phang}
{cmd:. use "didc_sim1.dta", clear}

{phang}
{cmd:. didc_bounds y, runvar(z) time(t) pre(0) post(1) id(id) c1(0(0.5)3) c2(0(0.5)3)}

{pstd}
The output gives a grid of identified sets. A lower bound that stays above zero
for a large part of the grid means the qualitative conclusion survives
substantial departures from time-invariant confounding.

{title:A note on the paper}

{pstd}
Lemma 6 is stated under "Assumptions 1, 2 and 5". The derivation needs
Assumption 3 (sharp discontinuities) and Assumption 6 (bounded variation), and
replaces Assumption 4. {cmd:didc_bounds} follows the derivation, not the printed
list. This is recorded in {cmd:docs/research-notes.md}, Section 9.

{title:References}

{pstd}
Manski, C. F., and J. V. Pepper (2018). How do right-to-carry laws affect crime
rates? Coping with ambiguous findings using partial identification.
{it:Review of Economics and Statistics} 100(2): 232-244.

{pstd}
Picchetti, P., C. C. X. Pinto, and S. T. Shinoki (2026).
Difference-in-discontinuities: estimation, inference and validity tests.
{it:arXiv preprint} arXiv:2405.18531.

{title:Also see}

{pmore}
Manual:  {help didc}, {help didc_test}
