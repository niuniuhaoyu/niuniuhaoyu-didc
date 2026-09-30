{smcl}
{* *! version 0.1.0  30sep2026}{...}
{vieweralsosee "didc" "help didc"}{...}
{title:Title}

{phang}
{bf:didc_test} --- Validity tests for the difference-in-discontinuities design

{title:Syntax}

{phang}
{cmd:didc_test} {it:depvar} {ifin} {cmd:,} {cmd:runvar(}{it:varname}{cmd:)}
{cmd:time(}{it:varname}{cmd:)} {cmdab:test(}{it:wald}{cmd:|}{it:ks}{cmd:)}
{cmd:pre(}{it:numlist}{cmd:)}
{hang2}{cmdab:cutoff(}{it:#}{cmd:)} {cmdab:p(}{it:#}{cmd:)}
{cmdab:kernel(}{it:kern}{cmd:)} {cmdab:bw(}{it:#}{cmd:)}
{cmdab:degree(}{it:#}{cmd:)} {cmdab:reps(}{it:#}{cmd:)} {cmdab:seed(}{it:#}{cmd:)}
{cmdab:level(}{it:#}{cmd:)} {cmd:noisily}

{title:Description}

{pstd}
{cmd:didc_test} implements the two validity tests of Section 4 of Picchetti,
Pinto & Shinoki (2026). Both ask whether the DiDC design is credible in the data
at hand by examining the pre-treatment periods, where no treatment of interest
was active.

{pstd}
{cmd:test(wald)} tests whether the discontinuity at the cutoff was the same in
every pre-treatment period. If a confounder is time-invariant, its jump does not
move; if it moves, the DiDC estimator does not identify the effect of interest.
The test is a joint Wald test on a stacked local-linear RD, and its bandwidth
defaults to the smallest period-specific CCT bandwidth, following the conclusion
of Appendix D.1 that any data-driven bandwidth appropriate to that RD works
well.

{pstd}
{cmd:test(ks)} tests whether the conditional mean of the outcome, {it:E[Y|Z]},
was the same in the two most recent pre-treatment periods. It uses a two-sample
Kolmogorov-Smirnov statistic on a polynomial basis, with a bootstrap that
imposes the null of a common conditional mean. This is the test of Section 4.2,
adapted from Whang (2001).

{pstd}
Rejection by either test is a warning about Assumption 4 (time-invariance of the
confounding effect), which is the assumption that the whole design rests on.
When it fails, {help didc_bounds} can report an identified set instead of a
point estimate.

{title:Options}

{phang}
{cmd:runvar(}{it:varname}{cmd:)}, {cmd:time(}{it:varname}{cmd:)},
{cmd:pre(}{it:numlist}{cmd:)} and {cmd:cutoff(}{it:#}{cmd:)} are as in
{help didc}.

{phang}
{cmd:test(wald|ks)} selects the test. Required.

{phang}
{cmd:pre(}{it:numlist}{cmd:)} must contain {bf:at least two} periods, all of
which must be pre-treatment. For {cmd:test(ks)} the two most recent of them are
compared; for {cmd:test(wald)} all of them enter the stacked regression.

{phang}
{cmd:p(}{it:#}{cmd:)} is the order of the local linear fit inside the stacked RD
(a hint) and {cmd:kernel(}{it:kern}{cmd:)} is its kernel. {cmd:bw(}{it:#}{cmd:)}
overrides the default bandwidth.

{phang}
{cmd:degree(}{it:#}{cmd:)} is the degree of the polynomial basis used by
{cmd:test(ks)}; default 3.

{phang}
{cmd:reps(}{it:#}{cmd:)} and {cmd:seed(}{it:#}{cmd:)} control the bootstrap of
{cmd:test(ks)}. Defaults are 999 and 20260930, so results are reproducible
without setting a seed.

{phang}
{cmd:noisily} shows the underlying regression for {cmd:test(wald)}.

{title:Stored results}

{pstd}
For {cmd:test(wald)}:

{p2colset 11 20 22 2}{...}
{p2col:{bf:r(F)}, {bf:r(df)}, {bf:r(p)}:joint test that the level jumps {it:gamma}{sub:k} are equal}
{p2col:{bf:r(F_theta)}, {bf:r(df_theta)}, {bf:r(p_theta)}:same for the slope changes {it:theta}{sub:k}}
{p2col:{bf:r(F_joint)}, {bf:r(df_joint)}, {bf:r(p_joint)}:both at once}
{p2col:{bf:r(gamma)}:1 x K matrix of level jumps, one column per pre period}
{p2col:{bf:r(theta)}:1 x K matrix of slope changes}
{p2col:{bf:r(bw)}, {bf:r(N)}:bandwidth and sample size used}
{p2col}

{pstd}
For {cmd:test(ks)}:

{p2colset 11 20 22 2}{...}
{p2col:{bf:r(ks_above)}, {bf:r(ks_below)}:KS statistics, above and below the cutoff}
{p2col:{bf:r(p_above)}, {bf:r(p_below)}:bootstrap p-values}
{p2col:{bf:r(N_above)}, {bf:r(N_below)}:observations on each side}
{p2col:{bf:r(period_new)}, {bf:r(period_old)}:the two periods compared}
{p2col}

{title:A note on equation (7) of the paper}

{pstd}
The paper prints the stacked regression of Section 4.1 as

{pmore}
{it:y} = sum over {it:k} of { [ {it:beta}{sub:-k}({it:Z} - {it:z0}) +
{it:theta}{sub:-k} {it:D}({it:Z} - {it:z0}) ] {it:T}{sub:-k} }

{pstd}
which has no level term in {it:D} and therefore only captures changes in the
{hi:slope} at the cutoff. A test of whether a discontinuity is stable needs the
level jump, so {cmd:didc_test} estimates

{pmore}
{it:y} = sum over {it:k} of { [ {it:alpha}{sub:k} +
{it:gamma}{sub:k} {it:D} + {it:beta}{sub:k}({it:Z} - {it:z0}) +
{it:theta}{sub:k} {it:D}({it:Z} - {it:z0}) ] {it:T}{sub:-k} }

{pstd}
and reports both. The headline test is on the level jumps {it:gamma}{sub:k},
because those are what "the discontinuity in period {it:k}" means;
{it:theta}{sub:k} is reported alongside. Simulations confirm that the test has
correct size under time-invariant confounding and near-unit power when the jump
moves, while the slope test has little power in that direction.

{title:A note on the variance of the Wald test}

{pstd}
The paper estimates the stacked regression by weighted least squares with
kernel weights and tests the restriction with the usual {it:F} statistic. If
those weights are supplied to Stata as {cmd:aweights}, Stata computes the
variance under the assumption that they are {hi:inverse-variance} weights, that
is, as if Var({it:u}{sub:i}) = {sigma}{sup:2}/{it:w}{sub:i}. The kernel weights
of a local polynomial regression are not inverse-variance weights: under
homoskedastic errors they are fixed weights on the estimator, and the correct
variance is {sigma}{sup:2}({it:X}'{it:WX}){sup:-1} with
{sigma}{sup:2} = sum({it:w}{sub:i}{it:u}{sub:i}{sup:2}) / (sum({it:w}{sub:i}) -
{it:k}). Using the {cmd:aweights} version understates the variance by roughly
the factor mean({it:w}), which in a local polynomial fit is well below one.

{pstd}
Measured on the paper's own model 1 design, with the default bandwidth, the
naive {cmd:aweights} test rejects a true null about 14% of the time at a
nominal 5%. {cmd:didc_test} therefore builds the variance from
{cmd:matrix accum} of {it:X}'{it:WX} and forms the Wald statistic itself, which
restores the nominal size. This is a deviation from the paper's implementation,
in the direction of the paper's stated test.

{title:Examples}

{phang}
{cmd:. use "didc_sim_multi.dta", clear}

{phang}
{cmd:. didc_test y, runvar(z) time(t) test(wald) pre(-1 0)}

{phang}
{cmd:. didc_test y, runvar(z) time(t) test(ks) pre(-1 0) reps(999)}

{title:References}

{pstd}
Hall, P., and J. L. Horowitz (1996). Bootstrap critical values for tests based on
generalized method of moments estimators. {it:Econometrica} 64(4): 891-916.

{pstd}
Picchetti, P., C. C. X. Pinto, and S. T. Shinoki (2026).
Difference-in-discontinuities: estimation, inference and validity tests.
{it:arXiv preprint} arXiv:2405.18531.

{pstd}
Whang, Y.-J. (2001). Consistent specification testing for conditional moment
restrictions. {it:Economics Letters} 71(3): 299-309.

{title:Also see}

{pmore}
Manual:  {help didc}, {help didc_bounds}
