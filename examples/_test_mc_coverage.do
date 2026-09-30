*==============================================================================
* _test_mc_coverage.do -- V4: Monte Carlo reproduction of Tables 1-4 of
* Picchetti, Pinto & Shinoki (2026) with tau = 0.
*
* For each of the four data-generating processes of Appendix D.2 the script
* reports average bias, RMSE and 95% coverage of the DiDC estimator, and, for
* model 1, the same quantities for a plain RD run on the post period alone.
* That last column is the point of the whole package: with a time-invariant
* confounder at the cutoff the plain RD is badly biased and its interval never
* covers, while didc is centred and its interval covers about 95% of the time.
*
*   do "D:/OpenCode/didc/examples/_test_mc_coverage.do" [reps] [nunits] [table]
*   defaults: reps 500, nunits 1000, table("reference/mc_results.csv")
*
* The paper reports 10,000 replications.  Raise reps to 10000 for that; the
* assertions are written to hold at either sample size.
*==============================================================================

clear all
set more off
adopath + "D:/OpenCode/didc"

global didc_npass = 0
global didc_nfail = 0

capture program drop _t
program define _t
    args label rcwant rcgot
    if `rcwant' == `rcgot' {
        display as result "  PASS  `label'   (rc = `rcgot')"
        global didc_npass = $didc_npass + 1
    }
    else {
        display as error "  FAIL  `label'   (expected rc `rcwant', got `rcgot')"
        global didc_nfail = $didc_nfail + 1
    }
end

*------------------------------------------------------------------ generator
* model 1: identical forms, time-invariant confounder c
* model 2: identical forms, no confounder
* model 3: time-varying forms, time-invariant confounder c
* model 4: time-varying forms, no confounder
capture program drop _mcgen
program define _mcgen
    args m n sd tau cjs
    clear
    set seed `sd'
    set obs `n'
    gen double id = _n
    gen double z  = 2*rbeta(2,4) - 1
    gen double beta_z = 0.48 + 1.27*z - 0.5*7.18*z^2 + 0.7*20.21*z^3 ///
                      + 1.1*21.54*z^4 + 1.5*7.33*z^5
    gen double R1 = 0.52 + 0.84*z - 0.1*3*z^2 - 0.3*7.99*z^3 ///
                  - 0.1*9.01*z^4 + 3.56*z^5
    gen double conf = cond(`m' == 2 | `m' == 4, 0, `cjs')
    gen double mu0 = .
    gen double mu1 = .
    replace mu0 = cond(z < 0, beta_z, R1 + conf)
    if `m' == 1 | `m' == 2 {
        replace mu1 = cond(z < 0, beta_z, R1 + conf + `tau')
    }
    else {
        replace mu1 = cond(z < 0, beta_z, 0.52 + 0.1*z + conf + `tau')
    }
    gen double y = .
    gen byte   t = .
    tempfile pre
    replace t = 0
    replace y = mu0 + rnormal(0, 0.1295)
    save `pre', replace
    replace t = 1
    replace y = mu1 + rnormal(0, 0.1295)
    append using `pre'
    sort id t
end

*------------------------------------------------------------------ arguments
local REPS  500
local NUNITS 1000
local TABLE "reference/mc_results.dta"
args a1 a2 a3
if "`a1'" != "" local REPS   `a1'
if "`a2'" != "" local NUNITS `a2'
if "`a3'" != "" local TABLE  "`a3'"

display as txt "Monte Carlo: `REPS' replications, n = `NUNITS', tau = 0"
display as txt "model |  av.bias      rmse    coverage   ci.length |  RD(post) bias  coverage"
display as txt "{hline 86}"

tempname pf
capture erase "`TABLE'"
postfile `pf' int model double bias double rmse double cover double cil ///
    double rdbias double rdcover using "`TABLE'", replace

forvalues m = 1/4 {

    local sum_t = 0
    local sum_t2 = 0
    local ncover = 0
    local sum_len = 0
    local sum_rd = 0
    local nrdcover = 0
    local nok = 0

    forvalues b = 1/`REPS' {
        quietly _mcgen `m' `NUNITS' `b' 0 0.5

        capture quietly didc y, runvar(z) time(t) pre(0) post(1) id(id) nolemma1
        if _rc == 0 {
            local ++nok
            local tt = e(tau_didc)
            local se = e(se_didc)
            local lo = e(ci_didc_l)
            local hi = e(ci_didc_r)
            local sum_t  = `sum_t'  + `tt'
            local sum_t2 = `sum_t2' + (`tt')^2
            local sum_len = `sum_len' + (`hi' - `lo')
            if `lo' <= 0 & 0 <= `hi' local ++ncover
        }

        * plain RD on the post period, the estimator the confounder breaks
        capture quietly rdrobust y z if t == 1, c(0)
        if _rc == 0 {
            local rd = e(tau_bc)
            local rlo = e(ci_l_rb)
            local rhi = e(ci_r_rb)
            local sum_rd = `sum_rd' + `rd'
            if `rlo' <= 0 & 0 <= `rhi' local ++nrdcover
        }
    }

    local bias    = `sum_t'/`nok'
    local rmse    = sqrt(`sum_t2'/`nok')
    local cover   = `ncover'/`nok'
    local cil     = `sum_len'/`nok'
    local rdbias  = `sum_rd'/`nok'
    local rdcover = `nrdcover'/`nok'

    display as txt "  `m'   | " as result %9.5f `bias' as txt %11.5f `rmse' ///
        "  " %8.4f `cover' "  " %8.4f `cil' as txt "  |  " ///
        as result %9.5f `rdbias' as txt "   " as result %8.4f `rdcover'

    post `pf' (`m') (`bias') (`rmse') (`cover') (`cil') (`rdbias') (`rdcover')

    if `m' == 1 {
        local ok = (abs(`bias') < 0.02)
        _t model1_bias_near_zero 1 `ok'
        local ok = (`cover' >= 0.90 & `cover' <= 0.98)
        _t model1_coverage_near_95 1 `ok'
        local ok = (abs(`rdbias') > 0.5)
        _t model1_plain_rd_is_badly_biased 1 `ok'
        local ok = (`rdcover' < 0.2)
        _t model1_plain_rd_does_not_cover 1 `ok'
    }
    if `m' == 2 {
        local ok = (abs(`bias') < 0.02)
        _t model2_bias_near_zero 1 `ok'
        local ok = (`cover' >= 0.88 & `cover' <= 0.99)
        _t model2_coverage_near_95 1 `ok'
    }
    if `m' == 3 {
        local ok = (abs(`bias') < 0.05)
        _t model3_bias_small_despite_time_varying_forms 1 `ok'
        local ok = (`cover' >= 0.88 & `cover' <= 0.99)
        _t model3_coverage_near_95 1 `ok'
    }
    if `m' == 4 {
        local ok = (abs(`bias') < 0.05)
        _t model4_bias_small 1 `ok'
        local ok = (`cover' >= 0.88 & `cover' <= 0.99)
        _t model4_coverage_near_95 1 `ok'
    }
}

postclose `pf'

display as txt "{hline 86}"
display as txt "  results written to `TABLE'"
display as txt "  NOTE: with a time-invariant confounder (model 1) the plain RD is biased"
display as txt "        by the size of the confounder and its interval never covers zero."

display as txt _n "{hline 60}"
display as txt "  passed: " as result "$didc_npass" as txt "    failed: " as result "$didc_nfail"
display as txt "{hline 60}"
if $didc_nfail > 0 {
    display as error "TEST SUITE FAILED"
    exit 9
}
display as result "ALL TESTS PASSED"
