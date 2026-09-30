*==============================================================================
* _test_engine.do -- v2: does the built-in Mata engine reproduce rdrobust?
*
* Two things are checked.
*
*   1. EXACTNESS of the point estimates.  At the same bandwidths, kernels and
*      polynomial orders, the built-in engine and rdrobust must agree on the
*      conventional estimate, the bias-corrected estimate, both side intercepts
*      and both side biases to machine precision.  This is an identity rather
*      than an approximation: both compute the same weighted least squares fit.
*
*   2. The VARIANCE conventions differ by design.  rdrobust reports CCT's
*      asymptotic variance with sigma^2 estimated by nearest neighbours; the
*      built-in engine reports the exact finite-sample variance of the same
*      linear functional.  The two are therefore compared by COVERAGE, in a
*      Monte Carlo, which is what a standard error is for.  Their ratio is also
*      printed, as information.
*
*   do "D:/OpenCode/didc/examples/_test_engine.do"
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

*==============================================================================
* Part 1: exactness of the point estimates
*==============================================================================
display as txt _n "== part 1: mata engine vs rdrobust at identical settings =="

* keep a copy of the panel; Delta Y is rebuilt by reshape, a different code
* path from the one inside _didc_prep
use "D:/OpenCode/didc/data/didc_sim1.dta", clear
keep id t z y
reshape wide y z, i(id) j(t)
gen double dy = y1 - y0
gen double zc = z1
keep dy zc
tempfile flat
save `flat', replace

local acc_ok = 1
local ncase  = 0
foreach case in 1 2 3 {

    if `case' == 1 local p 1
    if `case' == 1 local q 2
    if `case' == 1 local k "triangular"
    if `case' == 1 local h 0.4
    if `case' == 1 local b 0.6
    if `case' == 2 local p 1
    if `case' == 2 local q 2
    if `case' == 2 local k "uniform"
    if `case' == 2 local h 0.4
    if `case' == 2 local b 0.6
    if `case' == 3 local p 2
    if `case' == 3 local q 3
    if `case' == 3 local k "triangular"
    if `case' == 3 local h 0.5
    if `case' == 3 local b 0.8

    *---- rdrobust on the flat (dy, z) data ---------------------------
    use `flat', clear
    rename dy y
    rename zc z
    quietly rdrobust y z, c(0) p(`p') q(`q') kernel(`k') h(`h') b(`b')
    local rd_cl = e(tau_cl)
    local rd_bc = e(tau_bc)
    local rd_se_cl = e(se_tau_cl)
    local rd_se_rb = e(se_tau_rb)
    matrix BR = e(beta_Y_p_r)
    matrix BL = e(beta_Y_p_l)
    local rd_hp = BR[1,1]
    local rd_hm = BL[1,1]
    local rd_ba = e(bias_r)
    local rd_bb = e(bias_l)

    *---- the same thing through didc's built-in engine ----------------
    use "D:/OpenCode/didc/data/didc_sim1.dta", clear
    quietly didc y, runvar(z) time(t) pre(0) post(1) id(id) engine(mata) ///
        p(`p') q(`q') kernel(`k') h(`h') b(`b') nolemma1
    local my_cl = e(tau_didc)
    local my_bc = e(tau_didc_bc)
    local my_hp = e(mu_plus)
    local my_hm = e(mu_minus)
    local my_ba = e(bias_above)
    local my_bb = e(bias_below)
    local my_se = e(se_didc)

    local gap1 = max(abs(`rd_cl'-`my_cl'), abs(`rd_bc'-`my_bc'), ///
                     abs(`rd_hp'-`my_hp'), abs(`rd_hm'-`my_hm'))
    local gap2 = max(abs(`rd_ba'-`my_ba'), abs(`rd_bb'-`my_bb'))
    local ++ncase

    display as txt "  case `ncase': p=`p' q=`q' `k' h(`h') b(`b')"
    display as txt "     tau_cl  mata " %16.12f `my_cl' as txt "   rdrobust " %16.12f `rd_cl'
    display as txt "     tau_bc  mata " %16.12f `my_bc' as txt "   rdrobust " %16.12f `rd_bc'
    display as txt "     max |diff| over estimates and intercepts = " as result %12.4e `gap1'
    display as txt "     max |diff| over the two side biases     = " as result %12.4e `gap2'
    display as txt "     se_rb   mata " %16.12f `my_se' as txt "   rdrobust " %16.12f `rd_se_rb' ///
        as txt "   ratio " %6.4f (`my_se'/`rd_se_rb')

    if `gap1' > 1e-10 local acc_ok = 0
    if `gap2' > 1e-10 local acc_ok = 0
}
_t point_estimates_match_rdrobust 1 `acc_ok'

*-----------------------------------------------------------------------------
* p >= 2 is supported as well; both engines must still agree exactly
*-----------------------------------------------------------------------------
display as txt _n "== part 1b: p = 2 is accepted and still agrees exactly =="
use `flat', clear
rename dy y
rename zc z
quietly rdrobust y z, c(0) p(2) q(3) h(0.5) b(0.8)
local rd2 = e(tau_bc)
local rd2b = e(bias_l)

use "D:/OpenCode/didc/data/didc_sim1.dta", clear
capture quietly didc y, runvar(z) time(t) pre(0) post(1) id(id) ///
    engine(mata) p(2) q(3) h(0.5) b(0.8) nolemma1
local rc = _rc
_t engine_mata_accepts_p2 0 `rc'
if `rc' == 0 {
    local my2  = e(tau_didc_bc)
    local my2b = e(bias_below)
    display as txt "  tau_bc at p=2: mata " %16.12f `my2' as txt "   rdrobust " %16.12f `rd2'
    display as txt "  bias_l at p=2: mata " %16.12f `my2b' as txt "   rdrobust " %16.12f `rd2b'
    local ok = (abs(`my2' - `rd2') < 1e-10)
    _t p2_tau_bc_agrees 1 `ok'
    local ok = (abs(`my2b' - `rd2b') < 1e-10)
    _t p2_below_bias_agrees 1 `ok'
}

*==============================================================================
* Part 2: coverage of the built-in engine's interval
*==============================================================================
display as txt _n "== part 2: Monte Carlo coverage, model 1, tau = 0 =="

capture program drop _mce
program define _mce
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
    gen byte t = .
    tempfile pre
    replace t = 0
    replace y = mu0 + rnormal(0, 0.1295)
    save `pre', replace
    replace t = 1
    replace y = mu1 + rnormal(0, 0.1295)
    append using `pre'
    sort id t
end

local REPS 150
local cov_mata = 0
local cov_rd   = 0
local nmata = 0
local nrd   = 0

forvalues b = 1/`REPS' {
    quietly _mce 1 1000 `b' 0 0.5

    capture quietly didc y, runvar(z) time(t) pre(0) post(1) id(id) ///
        engine(mata) nolemma1
    if _rc == 0 {
        local ++nmata
        if e(ci_didc_l) <= 0 & 0 <= e(ci_didc_r) local ++cov_mata
    }

    capture quietly didc y, runvar(z) time(t) pre(0) post(1) id(id) nolemma1
    if _rc == 0 {
        local ++nrd
        if e(ci_didc_l) <= 0 & 0 <= e(ci_didc_r) local ++cov_rd
    }
}

local cm = `cov_mata'/`nmata'
local cr = `cov_rd'/`nrd'
display as txt "  `REPS' replications, n = 1000, tau = 0"
display as txt "    built-in engine (mata) : " as result %6.4f `cm' as txt "   (n = `nmata')"
display as txt "    rdrobust engine        : " as result %6.4f `cr' as txt "   (n = `nrd')"

local ok = (`cm' >= 0.86 & `cm' <= 1.0)
_t mata_engine_coverage_reasonable 1 `ok'
local ok = (`cr' >= 0.86 & `cr' <= 1.0)
_t rdrobust_engine_coverage_reasonable 1 `ok'

display as txt _n "{hline 60}"
display as txt "  passed: " as result "$didc_npass" as txt "    failed: " as result "$didc_nfail"
display as txt "{hline 60}"
if $didc_nfail > 0 {
    display as error "TEST SUITE FAILED"
    exit 9
}
display as result "ALL TESTS PASSED"
