*==============================================================================
* _test_validity_size.do -- V5/V6/V7: size and power of didc_test
*
* Data: three periods, t = -1 and t = 0 are pre-treatment, t = 1 is post.
*   NULL        the confounding jump is the same in both pre periods and the
*               functional form is the same  -> neither test should reject
*   ALT-LEVEL   the confounding jump differs between the pre periods
*               -> the Wald test should reject
*   ALT-SHAPE   the functional form above the cutoff differs between the pre
*               periods -> the KS test should reject
*
*   do "D:/OpenCode/didc/examples/_test_validity_size.do"
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
capture program drop _mk3
program define _mk3
    * args: c_m1 c_0 shape_m1 n seed
    args c_m1 c_0 shape_m1 n sd
    clear
    set seed `sd'
    set obs `n'
    gen double id = _n
    gen double z  = 2*rbeta(2,4)-1
    gen double beta_z = 0.48 + 1.27*z - 0.5*7.18*z^2 + 0.7*20.21*z^3 ///
                      + 1.1*21.54*z^4 + 1.5*7.33*z^5
    gen double R1 = 0.52 + 0.84*z - 0.1*3*z^2 - 0.3*7.99*z^3 ///
                  - 0.1*9.01*z^4 + 3.56*z^5
    gen double mu = .
    gen double y  = .
    gen byte   t  = .
    tempfile f1 f0
    replace t  = -1
    replace mu = cond(z<0, beta_z, cond(`shape_m1'==1, 0.52+0.1*z, R1) + `c_m1')
    replace y  = mu + rnormal(0,0.1295)
    save `f1', replace
    replace t  = 0
    replace mu = cond(z<0, beta_z, R1 + `c_0')
    replace y  = mu + rnormal(0,0.1295)
    save `f0', replace
    replace t  = 1
    replace mu = cond(z<0, beta_z, R1 + `c_0' + 0.3)
    replace y  = mu + rnormal(0,0.1295)
    append using `f0'
    append using `f1'
    sort id t
end

local NMCW = 200      // Monte Carlo replications for the Wald test
local NMCK = 100      // Monte Carlo replications for the KS test (slower)
local NUNITS = 1000

*-----------------------------------------------------------------------------
* V5: size of the Wald test under the null
*-----------------------------------------------------------------------------
display as txt _n "== V5: Wald test size under time-invariant confounding =="
local rej = 0
forvalues b = 1/`NMCW' {
    quietly _mk3 0.5 0.5 0 `NUNITS' `b'
    capture quietly didc_test y, runvar(z) time(t) test(wald) pre(-1 0)
    if _rc == 0 {
        if r(p) < 0.05 local ++rej
    }
}
local size = `rej'/`NMCW'
display as txt "    empirical size (nominal 5%) = " as result %6.4f `size'
local ok = (`size' >= 0.01 & `size' <= 0.11)
_t wald_size_within_tolerance 1 `ok'

*-----------------------------------------------------------------------------
* V6: power of the Wald test when the confounding jump changes
*-----------------------------------------------------------------------------
display as txt _n "== V6: Wald test power under time-varying confounding =="
local rej = 0
forvalues b = 1/`NMCW' {
    quietly _mk3 0.2 0.5 0 `NUNITS' `b'
    capture quietly didc_test y, runvar(z) time(t) test(wald) pre(-1 0)
    if _rc == 0 {
        if r(p) < 0.05 local ++rej
    }
}
local pw = `rej'/`NMCW'
display as txt "    power against a 0.3 change in the jump = " as result %6.4f `pw'
local ok = (`pw' > 0.5)
_t wald_power_above_half 1 `ok'

*-----------------------------------------------------------------------------
* V7: size and power of the KS test
*-----------------------------------------------------------------------------
display as txt _n "== V7: KS test size under a time-invariant conditional mean =="
local rej = 0
forvalues b = 1/`NMCK' {
    quietly _mk3 0.5 0.5 0 `NUNITS' `b'
    capture quietly didc_test y, runvar(z) time(t) test(ks) pre(-1 0) reps(199)
    if _rc == 0 {
        if r(p_above) < 0.05 local ++rej
    }
}
local sizesk = `rej'/`NMCK'
display as txt "    empirical size above the cutoff = " as result %6.4f `sizesk'
local ok = (`sizesk' <= 0.15)
_t ks_size_within_tolerance 1 `ok'

display as txt _n "== V7b: KS test power when the functional form changes =="
local rej = 0
forvalues b = 1/`NMCK' {
    quietly _mk3 0.5 0.5 1 `NUNITS' `b'
    capture quietly didc_test y, runvar(z) time(t) test(ks) pre(-1 0) reps(199)
    if _rc == 0 {
        if r(p_above) < 0.05 local ++rej
    }
}
local pwsk = `rej'/`NMCK'
display as txt "    power above the cutoff = " as result %6.4f `pwsk'
local ok = (`pwsk' > 0.5)
_t ks_power_above_half 1 `ok'

display as txt _n "{hline 60}"
display as txt "  passed: " as result "$didc_npass" as txt "    failed: " as result "$didc_nfail"
display as txt "{hline 60}"
if $didc_nfail > 0 {
    display as error "TEST SUITE FAILED"
    exit 9
}
display as result "ALL TESTS PASSED"
