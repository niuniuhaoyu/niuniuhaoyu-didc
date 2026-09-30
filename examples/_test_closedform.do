*==============================================================================
* _test_closedform.do -- V1: with Delta Y exactly linear on each side of the
* cutoff, the local linear estimator is exact and INDEPENDENT of the
* bandwidth.  This isolates the arithmetic from the bandwidth selector.
*
*   do "D:/OpenCode/didc/examples/_test_closedform.do"
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

*-----------------------------------------------------------------------------
* DGP: dy = 0.4 + 0.8*z + tau*(z>=0), exactly linear on each side, no noise.
* With p(1) any bandwidth reproduces tau exactly.
*-----------------------------------------------------------------------------
set seed 20260930
set obs 3000
gen double id = _n
gen double z  = 2*rbeta(2,4) - 1
local TAU = 0.3
gen double dy = 0.4 + 0.8*z + `TAU'*(z >= 0)
gen byte t = .
expand 2
bysort id: replace t = _n - 1
gen double y = .
replace y = dy - `TAU'*(z >= 0)         if t == 0
replace y = dy                            if t == 1
replace y = 0                             // placeholder, overwritten next line
replace y = dy*(t == 1) + (dy - `TAU'*(z >= 0))*(t == 0)
drop dy

display as txt _n "== V1: bandwidth invariance under exact linearity =="

* primary run with explicit bandwidths
capture noisily didc y, runvar(z) time(t) pre(0) post(1) id(id) ///
    h(0.2) b(0.3) p(1) q(2)
local rc = _rc
_t didc_runs 0 `rc'

if `rc' == 0 {
    local t1  = e(tau_didc)
    local t1b = e(tau_didc_bc)
    local mp1 = e(mu_plus)
    local mm1 = e(mu_minus)
    local l1  = e(lemma1_ok)
    display as txt "        tau_didc = " %12.8f `t1' "   tau_didc_bc = " %12.8f `t1b'
    display as txt "        mu_plus  = " %12.8f `mp1' "   mu_minus    = " %12.8f `mm1'
    display as txt "        lemma1_ok = `l1'"

    * exact recovery of tau
    local ok = (abs(`t1' - `TAU') < 1e-10)
    _t tau_recovered 1 `ok'
    local ok = (abs(`t1b' - `TAU') < 1e-10)
    _t taubc_recovered 1 `ok'
    local ok = (abs((`mp1' - `mm1') - `t1') < 1e-10)
    _t mu_plus_minus_eq_tau 1 `ok'
    local ok = (`l1' == 1)
    _t lemma1_check_pass 1 `ok'

    * now a different bandwidth: the answer must not move
    capture noisily didc y, runvar(z) time(t) pre(0) post(1) id(id) ///
        h(0.4) b(0.6) p(1) q(2)
    local rc2 = _rc
    _t didc_runs_again 0 `rc2'
    if `rc2' == 0 {
        local t2 = e(tau_didc)
        local gap = abs(`t1' - `t2')
        display as txt "        tau(h=0.4) = " %12.8f `t2' "   |gap| = " %12.10f `gap'
        local ok = (`gap' < 1e-10)
        _t bandwidth_invariance 1 `ok'
    }
}

*-----------------------------------------------------------------------------
* p(2) on a quadratic-in-z Delta Y: exact as well, checks the higher order path
*-----------------------------------------------------------------------------
display as txt _n "== V1b: p(2) on a quadratic differenced outcome =="
clear
set seed 7
set obs 3000
gen double id = _n
gen double z  = 2*rbeta(2,4) - 1
gen double dy = 0.4 + 0.8*z + 0.3*z^2 + 0.5*(z >= 0)
expand 2
bysort id: gen byte t = _n - 1
gen double y = .
replace y = dy                      if t == 1
replace y = dy - 0.5*(z >= 0)       if t == 0
drop dy

capture noisily didc y, runvar(z) time(t) pre(0) post(1) id(id) p(2) q(3) h(0.3) b(0.5)
local rc = _rc
_t didc_quad_runs 0 `rc'
if `rc' == 0 {
    local t = e(tau_didc)
    local d = abs(`t' - 0.5)
    display as txt "        tau_didc = " %12.8f `t' "  (true 0.5)   |gap| = " %14.12f `d'
    local ok = (`d' < 1e-8)
    _t quadratic_recovered 1 `ok'
}

display as txt _n "{hline 60}"
display as txt "  passed: " as result "$didc_npass" as txt "    failed: " as result "$didc_nfail"
display as txt "{hline 60}"
if $didc_nfail > 0 {
    display as error "TEST SUITE FAILED"
    exit 9
}
display as result "ALL TESTS PASSED"
