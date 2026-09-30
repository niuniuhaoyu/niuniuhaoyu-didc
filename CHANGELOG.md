# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.1] - 2026-09-30

Two errors in the built-in engine, both of which were **invisible at `p(1)`** —
which is exactly why the `1e-16` agreement with `rdrobust` failed to catch
them. The point estimates were correct throughout; only the variance and the
`p >= 2` bias correction were wrong, and the engine refused `p >= 2` in 0.2.0,
so no released number was affected.

### Fixed

- **The bias constant carried a spurious `factorial(p)`.** The leading bias is
  `h^(p+1) [e_0' Gamma^-1 vartheta_{p+1}] mu^(p+1)(0) / (p+1)!`, and because
  `mu^(p+1)(0) = (p+1)! beta_{p+1} / b^(p+1)` the two factorials cancel and no
  extra constant belongs there. `factorial(1) = 1`, so the `p(1)` validation
  could not see it; at `p(2)` it doubled the bias. Measured at `p(2), q(3)`,
  `h(0.5), b(0.8)`: the above-side bias was `0.004767` against `rdrobust`'s
  `0.002384`, and is now `0.00238360111601` against `0.00238360111601`.
- **The below-side bias carried a spurious orientation sign.** A Taylor
  expansion in the oriented coordinate `u = -(z-c)/h` suggests a factor
  `(-1)^(p+1)`, and that is what the engine applied. `(-1)^2 = 1`, so the `p(1)`
  validation could not see it either. Empirically `rdrobust` applies no
  orientation sign at `p(1)` or `p(2)`; the factor is removed. At `p(2)` the
  below-side bias is now `-0.001439927893` against `rdrobust`'s
  `-0.001439927893`.
- **The built-in engine's robust standard error was mis-scaled.** The weight
  vector of the bias-corrected estimator omitted the `1/b^(p+1)` factor that
  maps the `q`-order fit's leading coefficient onto the derivative at the
  cutoff. That factor cancels out of the point estimate — which is why the
  estimates were right and the mistake was invisible — but it does not cancel
  out of the variance. Without it the robust standard error was 0.83 times
  `rdrobust`'s and 95% coverage was 0.920.

### Changed

- **`p(2)` and above are now supported by the built-in engine**, at machine
  precision: with `p(2), q(3)` the two engines agree on the conventional
  estimate, the corrected estimate, both intercepts and both side biases to
  `8.3e-16`, `1.3e-15` for the biases alone.

| measurement, model 1, `n = 1000`, 150 replications | 0.2.0 | 0.2.1 |
|---|---|---|
| robust standard error relative to `rdrobust` | 0.83 | **0.94** |
| 95% coverage, built-in engine | 0.920 | **0.953** |
| 95% coverage, `rdrobust` engine | 0.960 | 0.960 |

### Added

- `examples/_test_engine.do` now covers `p(2)` as a case and prints the ratio
  of the two engines' robust standard errors and the built-in engine's coverage
  on every run. A mistake that leaves the point estimates exact can only be
  caught by coverage, so those two numbers are the guard.

## [0.2.0] - 2026-09-30

### Added

- **A built-in estimation engine.** `didc ..., engine(mata)` estimates the local
  polynomial fit, the bias correction and the variance without calling
  `rdrobust`, using a kernel written in Mata (`_didc_mata.ado`). With `h()` and
  `b()` supplied, no external package is involved at all; without them,
  `rdrobust` is used for bandwidth selection only. `engine(rdrobust)` remains
  the default and behaves exactly as in 0.1.0.
- `examples/_test_engine.do`, which compares the two engines quantity by
  quantity and measures their coverage by Monte Carlo. Added to the suite run
  by `examples/_test_all.do`.
- `e(bias_above)`, `e(bias_below)`, `e(engine)`, `e(kernel)`, `e(p_used)`,
  `e(q_used)`, `e(level_used)`.

### Design notes

- Every quantity the estimator needs is a linear functional of the outcome with
  known weights, so the built-in engine computes the variance directly from
  those weights, `Var(tau) = sum_s s2_s sum_i w_si^2`. That is the **exact**
  finite-sample variance for a fixed design and it avoids CCT's closed-form
  asymptotic algebra altogether.
- **Point estimates agree exactly.** At the same bandwidths, orders and kernel,
  the built-in engine reproduces `rdrobust`'s conventional estimate, corrected
  estimate, both side intercepts and both side biases to `2.8e-16` for the
  triangular kernel, `5.0e-16` for the uniform kernel and `2.8e-16` for
  epanechnikov, over six quantities at once.
- **The variance does not, and that is documented rather than hidden.**
  `rdrobust` reports CCT's asymptotic variance with `sigma^2` estimated by
  nearest neighbours; the built-in engine reports the exact finite-sample
  variance. Its interval is narrower, with a robust standard error about
  0.83 times `rdrobust`'s on the paper's model 1, and coverage of 0.920
  against 0.960 over 150 replications at `n = 1000`. That is why `rdrobust`
  stays the default. An attempt to back out the `sigma^2` implied by
  `rdrobust` from its conventional and robust standard errors produced a
  negative value, which shows that its robust variance is not a weighted sum
  of the two one-sided variances but CCT's closed form, including the extra
  term for estimating the bias.
  *Corrected in 0.2.1: the 0.83 ratio and the 0.920 coverage were caused by a
  missing `1/b^(p+1)` scale factor in the variance weights, not by the choice
  of variance convention. With that factor restored the ratio is 0.94 and
  coverage is 0.953.*
- **`p(2)` and above are refused by the built-in engine** rather than
  approximated. The bias correction reproduces `rdrobust` exactly at `p(1)`; at
  `p(2), q(3)` it differs by about `6e-04`, most likely because CCT's bias
  formula uses theoretical rather than empirical kernel moments. `engine(mata)`
  therefore exits with a pointer to `engine(rdrobust)`, so no user can
  accidentally publish a slightly different number.
- `docs/specs/2026-09-30-didc-v2-design.md` records what v2a delivered, the
  measurements above, and what v2b still owes: closing the coverage gap,
  lifting the `p >= 2` restriction, writing a bandwidth selector, and extending
  the built-in engine to `design(rcs)`.

## [0.1.0] - 2026-09-30

First release. Difference-in-discontinuities estimation, validity testing and
partial identification for Stata, following Picchetti, Pinto & Shinoki (2026),
arXiv:2405.18531.

### Added

- `didc`: local polynomial estimation of the treatment effect in a design where
  a confounder is assigned at the same cutoff. Supports panel data (the outcome
  difference is formed within unit and a single RD is run on it) and repeated
  cross-sections (one RD per period with a common bandwidth). Bias-corrected
  robust inference, and a built-in check of the paper's Lemma 1, that the RD of
  the differences and the difference of the RDs coincide when a common
  bandwidth is used.
- `didc_test`, `test(wald)`: stacked-RD Wald test that the discontinuity is the
  same in every pre-treatment period. Reports period-specific level jumps and
  slope changes, with joint and individual tests.
- `didc_test`, `test(ks)`: two-sample Kolmogorov-Smirnov test that the
  conditional mean of the outcome is time-invariant, with a bootstrap that
  imposes the null. Implemented in Mata.
- `didc_bounds`: partially identified sets under bounded variation, over a grid
  of the two sensitivity parameters, with breakdown values and explicit
  detection of empty identified sets.
- Simulated data for the paper's four Monte Carlo designs, plus a three-period
  dataset for the validity tests (`examples/didc_simdata.do`).
- English `README.md`, help files for all three commands, `didc.pkg`,
  `stata.toc`, AGPL-3.0 `LICENSE`.
- A test suite that runs without R: exact recovery and bandwidth invariance on a
  linear differenced outcome; removal of an explicitly simulated confounder;
  an independent hand-rolled weighted-least-squares estimator to check the
  local polynomial intercepts; a `reshape`-based reconstruction of the
  differenced outcome; Monte Carlo bias and coverage for the paper's models
  1-4; size and power of both validity tests; and closed-form checks of the
  bounds including the empty-set behaviour.

### Design notes

- Point estimation, bandwidth selection, bias correction and the robust
  variance are **delegated to `rdrobust`**, which is what the source paper
  itself does. The contributions of this package are the design layer, the two
  validity tests, the partial identification, and the Lemma 1 self-check.
  Self-contained estimation in Mata is planned for v2.
- **The Wald test uses a corrected variance.** The paper's stacked regression is
  a weighted least squares fit with kernel weights, tested with the usual `F`
  statistic. Passing kernel weights to Stata as `aweights` makes Stata treat
  them as inverse-variance weights, which understates the variance by roughly
  `mean(w)`. Measured on the paper's model 1 design at the default bandwidth the
  naive version rejects a true null 14.3% of the time at a nominal 5%.
  `didc_test` builds the correct homoskedastic WLS variance,
  `sigma2 (X'WX)^-1` with `sigma2 = sum(w u^2)/(sum(w) - k)`, from
  `matrix accum`, which brings the measured size to 3.3%. The test is now
  somewhat conservative rather than somewhat liberal, which is the safe
  direction, and the conservatism is explained by the local linear fit being
  misspecified against the strongly curved data-generating process.
- Measured behaviour, 200-300 Monte Carlo replications on the paper's designs:
  the KS test has size 0.050 and power 1.000; the corrected Wald test has size
  0.025 and power 0.990. The DiDC estimator has bias -0.004 and 95% coverage
  0.945 on model 1 and -0.008 / 0.950 on model 3, against the paper's
  -0.002 / 0.944 and -0.001 / 0.937. On the same data a plain post-treatment RD
  has bias 0.541 and coverage 0.000 in models 1 and 3, and 0.041 / 0.860 in
  models 2 and 4 against the paper's 0.043 / 0.861 for model 2. The paper does
  not report the magnitude of the confounder in model 1, so the size of the
  plain RD's bias depends on the value chosen here (0.5).
- Decisions taken where the source paper is internally inconsistent or silent
  are recorded in `docs/research-notes.md`, Section 9, and summarised in the
  README. The substantive ones: equation (7) is estimated with the level term
  the paper omits; `KS_-` is implemented on the below-cutoff support that the
  paper's own definition requires, not on the range printed in the formula;
  Lemma 6 is implemented from its derivation rather than its stated assumption
  list; and the empty identified set at `c1 = c2 = 0`, which the paper's Table 8
  reports, is treated as a real output rather than a typo.
- The condition for an empty identified set is `Delta Y^- != 0`, **not**
  `Delta Y^- > 0` as an earlier draft of the design notes claimed; the sign is
  irrelevant. Corrected in `docs/research-notes.md`, Section 11.4.
