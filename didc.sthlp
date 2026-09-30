{smcl}
{* *! version 0.1.0  30sep2026}{...}
{vieweralsosee "[R] rdrobust" "mansection R:rdrobust"}{...}
{vieweralsosee "didc_test" "help didc_test"}{...}
{vieweralsosee "didc_bounds" "help didc_bounds"}{...}
{title:Title}

{phang}
{bf:didc} --- Difference-in-discontinuities estimation

{title:Syntax}

{phang}
{cmd:didc} {it:depvar} {ifin} {cmd:,} {cmd:runvar(}{it:varname}{cmd:)}
{cmd:time(}{it:varname}{cmd:)} {cmd:pre(}{it:numlist}{cmd:)} {cmd:post(}{it:numlist}{cmd:)}
{hang2}{cmdab:design(}{it:panel}{cmd:|}{it:rcs}{cmd:)} {cmdab:id(}{it:varname}{cmd:)}
{cmdab:cutoff(}{it:#}{cmd:)}
{hang2}{cmdab:p(}{it:#}{cmd:)} {cmdab:q(}{it:#}{cmd:)} {cmdab:kernel(}{it:kern}{cmd:)}
{cmdab:bwselect(}{it:rule}{cmd:)} {cmdab:h(}{it:#}{cmd:)} {cmdab:b(}{it:#}{cmd:)}
{cmdab:level(}{it:#}{cmd:)}
{hang2}{cmdab:graph(}{it:name}{cmd:)} {cmd:nolemma1}
{hang2}{cmdab:engine(}{it:rdrobust}{cmd:|}{it:mata}{cmd:)}

{title:Description}

{pstd}
{cmd:didc} estimates the effect of a treatment that is assigned by a cutoff in
a running variable, in a setting where a {it:confounder} is assigned by the same
cutoff. This is the difference-in-discontinuities (DiDC) design. Under the
assumptions of Picchetti, Pinto & Shinoki (2026) stated in {help didc_bounds}
and below, the effect of interest is identified, while a standard regression
discontinuity that uses only the post-treatment cross-section is not.

{pstd}
The estimator is a local polynomial regression {hi:of the difference in
outcomes over time}, evaluated at the cutoff from each side. This is the
left-hand side of the paper's Lemma 1; the right-hand side is the difference
between the two period-specific discontinuities. The two are equal in the
population but not in finite samples unless a common bandwidth is used, and
{cmd:didc} always uses a common bandwidth so that they coincide.

{title:Options}

{phang}
{cmd:runvar(}{it:varname}{cmd:)} is the running variable. It must be
{hi:time-invariant} within unit; {cmd:didc} verifies this and reports the units
that violate it.

{phang}
{cmd:time(}{it:varname}{cmd:)} identifies the period, and
{cmd:pre(}{it:numlist}{cmd:)} and {cmd:post(}{it:numlist}{cmd:)} give one value
each. The data may contain only those two periods; any other non-missing value
is an error.

{phang}
{cmd:design(panel)} (the default) requires that the same units appear in both
periods, and requires {cmd:id()}. The outcome change is formed within unit and
one RD is run on it.

{phang}
{cmd:design(rcs)} is for repeated cross-sections, where the two periods are
independent samples, and does not allow {cmd:id()}. The outcome change cannot be
formed within unit, so one RD is run per period with a {hi:common} bandwidth and
the two are differenced. The common bandwidth is taken from the post period
unless {cmd:h()} and {cmd:b()} are supplied.

{phang}
{cmd:id(}{it:varname}{cmd:)} identifies the unit. Required with
{cmd:design(panel)}; not allowed with {cmd:design(rcs)}.

{phang}
{cmd:cutoff(}{it:#}{cmd:)} is the value of the running variable at which
treatment status changes. Default is 0. The running variable is centred
internally, so {cmd:h()} and {cmd:b()} are on the centred scale.

{phang}
{cmd:p(}{it:#}{cmd:)} and {cmd:q(}{it:#}{cmd:)} are the orders of the local
polynomial used for estimation and for the bias correction; {cmd:q()} must
exceed {cmd:p()}. Defaults are 1 and 2.

{phang}
{cmd:kernel(}{it:kern}{cmd:)}, {cmd:bwselect(}{it:rule}{cmd:)},
{cmd:h(}{it:#}{cmd:)} and {cmd:b(}{it:#}{cmd:)} are passed through to
{help rdrobust}. {cmd:h()} and {cmd:b()} accept one or two values, for below
and above the cutoff respectively.

{phang}
{cmd:level(}{it:#}{cmd:)} sets the confidence level; default is 95.

{phang}
{cmd:graph(}{it:name}{cmd:)} draws an RD plot of the differenced outcome and
stores it under {it:name}. Available with {cmd:design(panel)} only.

{phang}
{cmd:nolemma1} skips the Lemma 1 self-check. The check costs three extra
regressions and adds nothing to the estimate; it is on by default because it is
the cheap way to notice that something has gone wrong.

{phang}
{cmd:engine(rdrobust|mata)} selects the estimation engine. The default,
{cmd:engine(rdrobust)}, delegates the local polynomial fit, the bias
correction and the variance to {help rdrobust}, exactly as the source paper
does.

{pstd}
{cmd:engine(mata)} uses a self-contained kernel written in Mata. When
{cmd:h()} and {cmd:b()} are supplied it does not call {cmd:rdrobust} at all;
otherwise {cmd:rdrobust} is used for bandwidth selection only. Its point
estimates and bias corrections reproduce {cmd:rdrobust} to machine precision at
{cmd:p(1)} for the triangular, uniform and epanechnikov kernels.

{pstd}
The two engines differ in how the variance is computed. {cmd:rdrobust} reports
CCT's asymptotic variance with {it:sigma}{sup:2} estimated by nearest
neighbours; the built-in engine reports the exact finite-sample variance of the
same linear functional, computed from the weights of the estimator. Its interval
is therefore narrower and covers a little less often in finite samples: on the
paper's model 1 with {cmd:n = 1000} and 150 replications, coverage is 0.920 for
{cmd:engine(mata)} against 0.960 for {cmd:engine(rdrobust)}. That is why
{cmd:rdrobust} remains the default.

{warning:{cmd:engine(mata)} supports {cmd:p(1)} only.} At {cmd:p(2)} the bias
correction differs from {cmd:rdrobust} by about {cmd:6e-04}, so the built-in
engine refuses the request rather than returning a number that is close but not
identical. Use {cmd:engine(rdrobust)} for higher-order polynomials.

{pstd}
The engine is recorded in {cmd:e(engine)}.

{title:Assumptions}

{pstd}
Following the paper's numbering (Section 2):

{pmore}
1. {bf:Continuity.} All potential outcomes are continuous in {it:Z} at
{it:z0}. Note that this allows the confounder to jump: the confounder's jump is
handled by assumption 4, not ruled out here.

{pmore}
2. {bf:Random assignment at the cutoff.} Potential outcomes are independent of
the confounder and the treatment at {it:Z} = {it:z0}.

{pmore}
3. {bf:Sharp discontinuities.} Treatment take-up is 1 above the cutoff and 0
below, for both the confounder and the treatment of interest.

{pmore}
4. {bf:Time-invariance of the confounding effect.} The confounding effect at
the cutoff is the same before and after the treatment is introduced. This is
the assumption that makes DiDC work and it is the one to test:
see {help didc_test}.

{pmore}
5. {bf:No interaction.} The effect of the treatment of interest does not depend
on the confounder. Only needed to interpret the estimand as an effect that
generalises beyond the confounded units.

{pstd}
Sharply identified designs with imperfect compliance (fuzzy DiDC) are not
supported in this version.

{title:Stored results}

{pstd}
{cmd:didc} stores the usual {cmd:rdrobust} results, plus:

{p2colset 11 20 22 2}{...}
{p2col:{bf:e(tau_didc)}:conventional DiDC estimate}
{p2col:{bf:e(tau_didc_bc)}:bias-corrected estimate}
{p2col:{bf:e(se_didc)}:robust standard error of the bias-corrected estimate}
{p2col:{bf:e(ci_didc_l)}, {bf:e(ci_didc_r)}:robust confidence interval}
{p2col:{bf:e(mu_plus)}:{it:dY}{sup:+}, the intercept above the cutoff}
{p2col:{bf:e(mu_minus)}:{it:dY}{sup:-}, the intercept below the cutoff}
{p2col:{bf:e(lemma1_ok)}:1 if the two parameterisations agree with a common bandwidth}
{p2col:{bf:e(tau_delta_of_rds)}:the difference-of-RDs parameterisation}
{p2col:{bf:e(design)}:{cmd:panel} or {cmd:rcs}}
{p2col:{bf:e(cutoff)}:the cutoff used}
{p2col:{bf:e(n_units)}, {bf:e(n_below)}, {bf:e(n_above)}:sample sizes}
{p2col:{bf:e(bw_h_l)}, {bf:e(bw_h_r)}:estimation bandwidth, below and above}
{p2col:{bf:e(bw_b_l)}, {bf:e(bw_b_r)}:bias bandwidth, below and above}
{p2col:{bf:e(n_h_l)}, {bf:e(n_h_r)}:effective observations within the bandwidth}
{p2col}
{warning:Do not look for the bandwidths in {cmd:e(b)}.} As in {cmd:rdrobust},
{cmd:e(b)} holds the estimates, so the bandwidths are in {cmd:e(bw_*)}. This
avoids a collision between "bandwidth b" and "coefficient vector b".

{title:Things to watch for}

{pstd}
{bf:The treated group can be below the cutoff.} In Grembi, Nannicini &
Troiano (2016) treatment applies to municipalities with fewer than 5,000
inhabitants, i.e. {it:below} the cutoff. Nothing here assumes otherwise, but
when reading the output remember that "above the cutoff" and "treated" are not
synonyms.

{pstd}
{bf:Missing outcomes.} A unit with a missing outcome in either period cannot
contribute a difference and is dropped, with a note giving the count. A unit
observed in only one period is an error in {cmd:design(panel)}, not a silent
drop.

{pstd}
{bf:{cmd:test(wald)} needs an explicit membership test.} Internally,
{cmd:inlist(z, -1 0)} would be silently wrong because {cmd:inlist()} takes
comma-separated arguments; the command builds the indicator explicitly.

{title:Examples}

{phang}
{cmd:. use "didc_sim1.dta", clear}

{phang}
{cmd:. didc y, runvar(z) time(t) pre(0) post(1) id(id)}

{phang}
{cmd:. didc y, runvar(z) time(t) pre(0) post(1) id(id) p(2) q(3) kernel(uniform)}

{phang}
{cmd:. didc y, runvar(z) time(t) pre(0) post(1) id(id) graph(doseplot)}

{title:References}

{pstd}
Calonico, S., M. D. Cattaneo, and R. Titiunik (2014). Robust nonparametric
confidence intervals for regression-discontinuity designs.
{it:Econometrica} 82(6): 2295-2326.

{pstd}
Grembi, V., T. Nannicini, and U. Troiano (2016). Do fiscal rules matter?
{it:American Economic Journal: Applied Economics} 8(3): 1-30.

{pstd}
Picchetti, P., C. C. X. Pinto, and S. T. Shinoki (2026).
Difference-in-discontinuities: estimation, inference and validity tests.
{it:arXiv preprint} arXiv:2405.18531.

{title:Also see}

{pmore}
Manual:  {help didc_test}, {help didc_bounds}, {help rdrobust}
