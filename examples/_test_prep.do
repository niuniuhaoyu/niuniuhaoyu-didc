*==============================================================================
* _test_prep.do -- tests for _didc_prep (validation and Delta Y construction)
*
*   do "D:/OpenCode/didc/examples/_test_prep.do"
*
* Labels must contain no spaces (args splits on spaces).
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

*---------------------------------------------------------------- case 5 ------
display as txt _n "== case 5: Delta Y direction and build correctness =="
clear
input byte id byte t double z double y
1 0  0.5  10
1 1  0.5  13
2 0 -0.5  20
2 1 -0.5  18
end
preserve
capture noisily _didc_prep build, depvar(y) runvar(z) time(t) pre(0) post(1) id(id)
local rc = _rc
local ndrop = r(n_dropped_missing_dy)
_t build_clean_panel 0 `rc'
if `rc' == 0 {
    local nrem = _N
    _t two_units_remain 2 `nrem'
    local ok = 1
    quietly count if id == 1 & _didc_dy ==  3
    if r(N) != 1 local ok = 0
    quietly count if id == 2 & _didc_dy == -2
    if r(N) != 1 local ok = 0
    _t dy_is_post_minus_pre 1 `ok'
    local ok = 1
    quietly count if _didc_zc == 0.5 & id == 1
    if r(N) != 1 local ok = 0
    _t zc_is_z_minus_cutoff 1 `ok'
    _t no_units_dropped 0 `ndrop'
}
restore

*---------------------------------------------------------------- case 1 ------
display as txt _n "== case 1: clean two-period panel passes check =="
use "D:/OpenCode/didc/data/didc_sim1.dta", clear
capture noisily _didc_prep check, depvar(y) runvar(z) time(t) pre(0) post(1) id(id)
local rc     = _rc
local nunits = r(n_units)
local nbelow = r(n_below)
local nabove = r(n_above)
_t check_didc_sim1 0 `rc'
if `rc' == 0 {
    display as txt "        n_units = `nunits'   below = `nbelow'   above = `nabove'"
    _t n_units_eq_2000 2000 `nunits'
}

*---------------------------------------------------------------- case 2 ------
display as txt _n "== case 2: three periods must be rejected =="
clear
input byte id byte t double z double y
1 0  0.5 1
1 1  0.5 2
1 2  0.5 3
2 0 -0.5 1
2 1 -0.5 2
2 2 -0.5 3
end
capture _didc_prep check, depvar(y) runvar(z) time(t) pre(0) post(1) id(id)
local rc = _rc
_t three_periods_rejected 198 `rc'

*---------------------------------------------------------------- case 3 ------
display as txt _n "== case 3: runvar varying over time must be rejected =="
clear
input byte id byte t double z double y
1 0  0.5 1
1 1  0.6 2
2 0 -0.5 1
2 1 -0.5 2
end
capture _didc_prep check, depvar(y) runvar(z) time(t) pre(0) post(1) id(id)
local rc = _rc
_t time_varying_runvar_rejected 459 `rc'

*---------------------------------------------------------------- case 4 ------
display as txt _n "== case 4: pre() == post() must be rejected =="
use "D:/OpenCode/didc/data/didc_sim1.dta", clear
capture _didc_prep check, depvar(y) runvar(z) time(t) pre(1) post(1) id(id)
local rc = _rc
_t pre_equals_post_rejected 198 `rc'

*---------------------------------------------------------------- case 6 ------
display as txt _n "== case 6: design()/id() consistency =="
use "D:/OpenCode/didc/data/didc_sim1.dta", clear
capture _didc_prep check, depvar(y) runvar(z) time(t) pre(0) post(1) design(panel)
local rc = _rc
_t panel_without_id_rejected 198 `rc'
capture _didc_prep check, depvar(y) runvar(z) time(t) pre(0) post(1) design(rcs) id(id)
local rc = _rc
_t rcs_with_id_rejected 198 `rc'
capture _didc_prep check, depvar(y) runvar(z) time(t) pre(0) post(1) design(bogus) id(id)
local rc = _rc
_t unknown_design_rejected 198 `rc'
capture noisily _didc_prep check, depvar(y) runvar(z) time(t) pre(0) post(1) design(rcs)
local rc = _rc
_t rcs_without_id_accepted 0 `rc'

*---------------------------------------------------------------- case 7 ------
display as txt _n "== case 7: units missing a period must be rejected =="
clear
input byte id byte t double z double y
1 0  0.5 1
1 1  0.5 2
2 0 -0.5 1
3 1 -0.3 2
end
capture _didc_prep check, depvar(y) runvar(z) time(t) pre(0) post(1) id(id)
local rc = _rc
_t one_period_units_rejected 459 `rc'

*---------------------------------------------------------------- case 8 ------
display as txt _n "== case 8: missing outcome is dropped, not fatal =="
clear
input byte id byte t double z double y
1 0  0.5 1
1 1  0.5 2
2 0 -0.5 1
2 1 -0.5 .
end
preserve
capture noisily _didc_prep build, depvar(y) runvar(z) time(t) pre(0) post(1) id(id)
local rc    = _rc
local nrem  = _N
local ndrop = r(n_dropped_missing_dy)
_t build_tolerates_missing_outcome 0 `rc'
if `rc' == 0 {
    _t one_unit_remains 1 `nrem'
    _t one_unit_dropped 1 `ndrop'
}
restore

*---------------------------------------------------------------- summary ----
display as txt _n "{hline 60}"
display as txt "  passed: " as result "$didc_npass" as txt "    failed: " as result "$didc_nfail"
display as txt "{hline 60}"
if $didc_nfail > 0 {
    display as error "TEST SUITE FAILED"
    exit 9
}
display as result "ALL TESTS PASSED"
