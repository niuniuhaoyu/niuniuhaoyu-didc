# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
