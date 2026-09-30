*==============================================================================
* _test_vs_wls.do -- V3: independent numerical check against a hand-rolled
* local linear estimator, plus an independent reconstruction of Delta Y.
*
* R is not installed on this machine, so instead of comparing with the R
* package rdrobust we re-derive the CONVENTIONAL local polynomial estimate from
* first principles:
*
*   on each side of the cutoff, fit  dy = a + b*z  by weighted least squares
*   with weights equal to the kernel, restricted to |z| <= h.  The intercept
*   a is mu_minus / mu_plus.  This uses regress, not rdrobust.
*
* Delta Y itself is rebuilt through reshape, i.e. a different code path from
* the one inside _didc_prep, so the plumbing is checked as well.
*
* When R is available, examples/reference/run_rdrobust.R reproduces the same
* comparison against the R package rdrobust (see README).
*
*   do "D:/OpenCode/didc/examples/_test_vs_wls.do"
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

display as txt _n "== V3: didc vs a hand-rolled local linear estimator =="

*---- didc ---------------------------------------------------------------
use "D:/OpenCode/didc/data/didc_sim1.dta", clear
capture noisily didc y, runvar(z) time(t) pre(0) post(1) id(id) p(1) q(2)
local rc = _rc
_t didc_runs 0 `rc'
if `rc' != 0 exit 9

local H        = e(bw_h_l)
local mu_plus  = e(mu_plus)
local mu_minus = e(mu_minus)
local tau      = e(tau_didc)
display as txt "    didc used h = " %8.5f `H' "   (mserd picks h_l = h_r = h)"

*---- independent reconstruction ----------------------------------------
* reshape wide is a completely different route to Delta Y
use "D:/OpenCode/didc/data/didc_sim1.dta", clear
keep id t z y
reshape wide y z, i(id) j(t)
quietly count if z0 != z1
_t runvar_time_invariant_on_reconstruction 0 r(N)
gen double dy_r = y1 - y0
gen double zc_r = z1

quietly count
local Nunits = r(N)
display as txt "    reconstructed N = `Nunits'"

*---- hand-rolled local linear, triangular kernel ------------------------
* weights: K(u) = (1 - |u|)+ with u = z/h ; only |z| <= h contributes
gen double w = max(0, 1 - abs(zc_r/`H')) if abs(zc_r) <= `H'

quietly regress dy_r zc_r if zc_r >= 0 [aw = w]
local mu_plus_w  = _b[_cons]
quietly regress dy_r zc_r if zc_r <  0 [aw = w]
local mu_minus_w = _b[_cons]
local tau_w      = `mu_plus_w' - `mu_minus_w'

display as txt "                          didc        hand-rolled WLS      |gap|"
display as txt "    mu_plus   " %14.8f `mu_plus'  %18.8f `mu_plus_w'  %14.10f abs(`mu_plus'-`mu_plus_w')
display as txt "    mu_minus  " %14.8f `mu_minus' %18.8f `mu_minus_w' %14.10f abs(`mu_minus'-`mu_minus_w')
display as txt "    tau_didc  " %14.8f `tau'      %18.8f `tau_w'      %14.10f abs(`tau'-`tau_w')

local ok = (abs(`mu_plus'  - `mu_plus_w')  < 1e-6)
_t mu_plus_matches_wls 1 `ok'
local ok = (abs(`mu_minus' - `mu_minus_w') < 1e-6)
_t mu_minus_matches_wls 1 `ok'
local ok = (abs(`tau'      - `tau_w')      < 1e-6)
_t tau_matches_wls 1 `ok'

*-----------------------------------------------------------------------------
* same exercise with a uniform kernel, to rule out a kernel-convention artefact
*-----------------------------------------------------------------------------
display as txt _n "== V3b: same check with kernel(uniform) =="
use "D:/OpenCode/didc/data/didc_sim1.dta", clear
capture noisily didc y, runvar(z) time(t) pre(0) post(1) id(id) ///
    p(1) q(2) kernel(uniform)
local rc = _rc
_t didc_uniform_runs 0 `rc'
if `rc' == 0 {
    local H2       = e(bw_h_l)
    local mu_plus2 = e(mu_plus)
    local mu_min2  = e(mu_minus)

    use "D:/OpenCode/didc/data/didc_sim1.dta", clear
    keep id t z y
    reshape wide y z, i(id) j(t)
    gen double dy_r = y1 - y0
    gen double zc_r = z1
    gen byte inside = (abs(zc_r) <= `H2')
    quietly regress dy_r zc_r if zc_r >= 0 & inside
    local mp2 = _b[_cons]
    quietly regress dy_r zc_r if zc_r <  0 & inside
    local mm2 = _b[_cons]

    display as txt "    mu_plus : didc " %12.8f `mu_plus2' "   WLS " %12.8f `mp2' ///
        "   |gap| " %12.10f abs(`mu_plus2'-`mp2')
    display as txt "    mu_minus: didc " %12.8f `mu_min2'  "   WLS " %12.8f `mm2' ///
        "   |gap| " %12.10f abs(`mu_min2'-`mm2')
    local ok = (abs(`mu_plus2'-`mp2') < 1e-6)
    _t uniform_mu_plus_matches_wls 1 `ok'
    local ok = (abs(`mu_min2'-`mm2') < 1e-6)
    _t uniform_mu_minus_matches_wls 1 `ok'
}

display as txt _n "{hline 60}"
display as txt "  passed: " as result "$didc_npass" as txt "    failed: " as result "$didc_nfail"
display as txt "{hline 60}"
if $didc_nfail > 0 {
    display as error "TEST SUITE FAILED"
    exit 9
}
display as result "ALL TESTS PASSED"
