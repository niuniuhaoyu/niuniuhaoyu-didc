# didc

**Difference-in-Discontinuities for Stata**

[![Stata 17+](https://img.shields.io/badge/Stata-17%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 0.1.0](https://img.shields.io/badge/Version-0.1.0-green.svg)](CHANGELOG.md)
[![Requires: rdrobust](https://img.shields.io/badge/requires-rdrobust-informational.svg)](https://rdpackages.github.io/rdrobust/)

`didc` implements **difference-in-discontinuities (DiDC)** estimation, inference,
validity testing and partial identification for Stata, following
Picchetti, Pinto & Shinoki (2026), *Difference-in-Discontinuities: Estimation,
Inference and Validity Tests* ([arXiv:2405.18531](https://arxiv.org/abs/2405.18531)).

---

## The problem it solves

In many policy settings the treatment of interest is assigned by a cutoff in a
running variable **and so is something else**. Italian municipalities below
5,000 inhabitants had their fiscal rules relaxed in 2001 — and also had a
higher mayoral wage, assigned by the same 5,000-inhabitant cutoff. Other
examples: geographic boundaries that host several programmes at once, school
thresholds that trigger both a class-size rule and a funding rule.

Standard regression discontinuity is invalid here: the running variable is not
continuous in mean potential outcomes at the cutoff, because the *confounder*
jumps there too. Difference-in-differences is usually not available either,
since the units on either side of the cutoff are not comparable in levels.

DiDC fixes this by using the **change over time**: a confounder that is
constant over time contributes the same jump in the pre and post periods and
cancels out of the difference.

## Installation

```stata
* dependency
ssc install rdrobust, replace

* this package
net install didc, from("https://raw.githubusercontent.com/niuniuhaoyu/niuniuhaoyu-didc/main/") replace
```

`didc` needs a clear rejection of a common failure mode: if `rdrobust` is
missing it says so and tells you how to install it, rather than failing with
`r(111)`.

## Quick start

```stata
use "didc_sim1.dta", clear          // simulated data shipped with the package
didc y, runvar(z) time(t) pre(0) post(1) id(id)

didc_test y, runvar(z) time(t) pre(-1 0) test(wald)   // is the confounder stable?
didc_test y, runvar(z) time(t) pre(-1 0) test(ks)     // is E[Y|Z] stable?

didc_bounds y, runvar(z) time(t) pre(0) post(1) id(id) c1(0(0.5)5) c2(0(0.5)5)
```

A complete, reproducible walkthrough is in `examples/didc_example.do`.

## Data layout

The command takes **long** data: one row per unit and period.

| `design()` | meaning | `id()` | what is estimated |
|---|---|---|---|
| `panel` (default) | the same units appear in both periods | required | `DeltaY = y(post) - y(pre)` is formed within unit, then **one** RD is run on `(DeltaY, Z)` |
| `rcs` | repeated cross-sections; the two periods are independent samples | not allowed | **two** RDs are run, one per period, with a common bandwidth, and differenced |

`runvar()` must be **time-invariant** within unit; the command checks this and
tells you which units violate it.

## Method

### Estimand

Let `Z` be the running variable, `z0 = cutoff()`, `D_{i,0}` the confounding
policy and `D_{i,1}` the treatment of interest, both assigned by `1{Z >= z0}`.
The target parameter is

```
tau_c = E[ Y_{i,1}(1,1) - Y_{i,1}(1,0) | Z_i = z0 ]
```

the effect of the treatment of interest at the cutoff for units exposed to the
confounder. Under Assumptions 1-4 of the paper it equals the DiDC estimand

```
tau_didc = tau_1^RD - tau_0^RD
```

### Lemma 1: the estimator is a single RD on the difference

The paper's key computational result is that the *difference of the RDs* equals
the *RD of the differences*:

```
tau_didc = lim_{z->z0+} E[ dY | Z=z ] - lim_{z->z0-} E[ dY | Z=z ],   dY = y_1 - y_0
```

so the estimator is a local polynomial fit of `dY` on `Z`, exactly as in
Cattaneo, Calonico & Titiunik (2014), which is what `rdrobust` does.

**A finite-sample subtlety that `didc` handles for you.** The identity is a
statement about population limits. If you estimate the two sides separately and
let each regression pick its own MSE-optimal bandwidth, the two
parameterisations *do not* agree numerically. On the simulated data:

| | bandwidth | estimate |
|---|---|---|
| RD of `dY` | 0.2165 | 0.300031 |
| `RD(y_post) - RD(y_pre)` | each its own | 0.290737 |
| both, common bandwidth | 0.2165 | 0.300031 vs 0.300031 |

`didc` therefore always reports both parameterisations **with a common
bandwidth**, and returns `e(lemma1_ok)` as a built-in check. No other Stata
package does this.

### Inference

Local polynomial estimation with the MSE-optimal plug-in bandwidth selector,
bias correction, and the robust confidence interval of Calonico, Cattaneo &
Titiunik (2014). All of that comes from `rdrobust`; see *Implementation notes*
below.

### Validity tests (`didc_test`)

| test | question | specification |
|---|---|---|
| `test(wald)` | Is the discontinuity the same in every pre-treatment period? | a stacked local-linear RD over the periods in `pre()`, with period-specific level jumps `gamma_k` and slope changes `theta_k`; joint Wald tests |
| `test(ks)` | Is `E[Y\|Z]` the same in the two most recent pre-treatment periods? | a two-sample Kolmogorov-Smirnov statistic on a polynomial basis, with a bootstrap that **imposes the null** |

Rejection of either undermines Assumption 4 (time-invariance of the confounding
effect) and therefore the DiDC design itself. Neither test exists elsewhere in
the Stata ecosystem.

**A variance fix inside `test(wald)` that matters.** The paper estimates the
stacked regression by weighted least squares with kernel weights and tests with
the usual `F` statistic. Handing kernel weights to Stata as `aweights` makes
Stata treat them as inverse-variance weights, which understates the variance by
roughly `mean(w)` and makes the test far too liberal — measured at about **14%**
rejection at a nominal 5% on the paper's model 1. `didc_test` therefore forms
the correct homoskedastic WLS variance,
`sigma2 (X'WX)^{-1}` with `sigma2 = sum(w u^2)/(sum(w) - k)`, from `matrix accum`
and builds the Wald statistic directly. This restores the nominal size and is a
deviation from the paper's implementation in the direction of the paper's stated
test.

### Partial identification (`didc_bounds`)

When Assumption 4 is not credible, replace it with bounded variation
(Assumption 6) and report an identified set instead of a point:

```
tau_c in [ max{ mu_plus - c1 , tau_didc - c2 } ,
           min{ mu_plus + c1 , tau_didc + c2 } ]
```

with `mu_plus = dY^+`. `c1` bounds a time trend in the confounded outcome;
`c2` bounds drift in the confounding effect itself. The command reports the
whole `c1 x c2` grid and the **breakdown values**: how large `c1` or `c2` can
get before you can no longer sign the effect.

**Empty identified sets are a feature, not a bug.** At `c1 = c2 = 0` the set is
empty whenever `dY^- != 0`, i.e. almost always in real data. That is direct
evidence against Assumption 4, and `didc_bounds` says so instead of erroring.

## Implementation notes (what is and is not original)

| component | how it is produced | original? |
|---|---|---|
| panel/rcs handling, `dY` construction, time-invariance check | written here | yes |
| point estimate, MSE-optimal bandwidth, bias correction, robust CI | **delegated to `rdrobust`** | no |
| Lemma 1 self-check | written here | yes |
| stacked-RD Wald test | written here | yes, no precedent in Stata |
| bootstrap KS test | written here | yes, no precedent in Stata |
| partial identification, breakdown values, empty-set detection | written here | yes, no precedent in Stata |

Delegating the estimation kernel is what the source paper does — it explicitly
builds on Calonico, Cattaneo & Titiunik (2014). Hand-rolling the bandwidth
selector and the robust variance would introduce differences that nobody could
verify, and would destroy the one thing this package can offer that a quick port
cannot: an **independent numerical check** of the estimates. See
`docs/specs/2026-09-30-didc-design.md`, section 5.1.

Self-contained estimation (bandwidth selector, local polynomial, bias
correction and robust variance implemented in Mata) is on the roadmap for v2;
when it lands, the two implementations check each other.

### Difference from `RDDID`

`RDDID` (SSC `s459586`, Jonathan Dries) performs a *cross-sectional* DiDC: the
discontinuity in the treated group minus the discontinuity in the control
group. It has no time dimension, no validity tests and no partial
identification. `didc` estimates the *time-differenced* design of Grembi et al.
(2016), which is a different research design. The two share `rdrobust` as a
computational backend and nothing else.

### Inconsistencies in the source paper and what `didc` does about them

Recorded in full in `docs/research-notes.md`, Section 9. The substantive ones:

1. **Equation (7)** is printed without the level term `gamma_k D`, which would
   leave only slope changes at the cutoff and no test of the discontinuity
   itself. `didc_test` estimates both and reports both; the headline test is on
   the level jumps.
2. **`KS_-`** is printed summing over `Z_i >= z0`, contradicting its own
   definition of the below-cutoff support. `didc_test` uses `Z_i < z0`,
   symmetric with `KS_+`.
3. **Lemma 6** is stated under "Assumptions 1, 2 and 5"; the derivation needs
   Assumption 3 and Assumption 6. `didc_bounds` follows the derivation.
4. **Tables 1-4** are headed "10,000 Monte Carlo experiments" while Section 6
   says 1,000. `examples/_test_mc_coverage.do` takes the replication count as
   an argument and works at either.

## Stored results

`didc` returns the usual `rdrobust` results plus:

| `e()` | contents |
|---|---|
| `e(tau_didc)` | conventional DiDC estimate |
| `e(tau_didc_bc)` | bias-corrected estimate |
| `e(se_didc)`, `e(ci_didc_l)`, `e(ci_didc_r)` | robust standard error and interval |
| `e(mu_plus)`, `e(mu_minus)` | `dY^+` and `dY^-`, the left/right intercepts |
| `e(lemma1_ok)` | 1 if both parameterisations agree with a common bandwidth |
| `e(tau_delta_of_rds)` | the difference-of-RDs parameterisation |
| `e(design)`, `e(cutoff)`, `e(n_units)`, `e(n_below)`, `e(n_above)` | design and sample |
| `e(bw_h_l)`, `e(bw_h_r)`, `e(bw_b_l)`, `e(bw_b_r)` | bandwidths (note: **not** `e(b)`, which holds the estimate) |

`didc_test` returns `r(p)`, `r(F)`, `r(df)`, the `gamma` and `theta` matrices
and the slope/joint tests for `test(wald)`, and `r(ks_above)`, `r(ks_below)`,
`r(p_above)`, `r(p_below)` for `test(ks)`.

`didc_bounds` returns `e(bounds_lb)`, `e(bounds_ub)`, `e(bounds_empty)`,
`e(breakdown_c1)`, `e(breakdown_c2)`.

## Evidence that the numbers are right

| check | file | result |
|---|---|---|
| exact recovery and bandwidth invariance on a linear differenced outcome | `examples/_test_closedform.do` | exact to `1e-10` |
| the confounder is removed: `didc` recovers `tau` while a post-period-only RD returns `c + tau` | `examples/_test_confounder.do` | as predicted |
| **independent** local-linear estimator (weighted least squares via `regress`, no `rdrobust`) | `examples/_test_vs_wls.do` | agrees to `0.0000000000` for both kernel types |
| rebuilt `dY` through `reshape` rather than the package's own path | `examples/_test_vs_wls.do` | identical |
| Monte Carlo bias and coverage, models 1-4 | `examples/_test_mc_coverage.do` | model 1: bias `-0.002`, coverage `0.944`, plain post-period RD bias `0.542` and coverage `0.000`; model 2: coverage `0.944`, plain RD `0.042` / `0.846` against the paper's `0.043` / `0.861` |
| size and power of both validity tests | `examples/_test_validity_size.do` | KS: size `0.050`, power `1.000`; Wald: size `0.033`, power `1.000` |
| closed-form bounds, empty sets, breakdown values | `examples/_test_bounds.do` | cell-by-cell exact |

Run the whole suite at once with:

```stata
do "examples/_test_all.do"        // or:  do "examples/_test_all.do" 1000
```

`examples/reference/run_rdrobust.R` reproduces the `_test_vs_wls.do` comparison
against the **R** package `rdrobust` on a machine where R is available; the
in-repository check does not need R.

## Citation

BibTeX for the method:

```bibtex
@article{picchetti2026didc,
  author  = {Picchetti, Pedro and Pinto, Cristine C. X. and Shinoki, St{\'e}phanie T.},
  title   = {Difference-in-Discontinuities: Estimation, Inference and Validity Tests},
  journal = {arXiv preprint arXiv:2405.18531},
  year    = {2026}
}
```

The empirical example follows:

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

## License

AGPL-3.0. See `LICENSE`.
