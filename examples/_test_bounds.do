*==============================================================================
* _test_bounds.do -- tests for didc_bounds, and an end-to-end check of
* design(rcs) against design(panel)
*
* The identified set is a closed form, so it is checked cell by cell:
*   LB(c1,c2) = max{ mu_plus - c1 , tau_didc - c2 }
*   UB(c1,c2) = min{ mu_plus + c1 , tau_didc + c2 }
*   EMPTY     = (LB > UB), which at c1 = c2 = 0 holds iff Delta Y^- != 0
* The assertions on the closed form use the ESTIMATED mu_plus and tau_didc, so
* they are exact identities rather than sampling statements.
*
*   do "D:/OpenCode/didc/examples/_test_bounds.do"
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
capture program drop _mkb
program define _mkb
    * args: level_below slope_above intercept_above jump n
    * dy is the intended per-unit DIFFERENCE y(post) - y(pre).  Levels are
    * built by putting pure noise in the pre period, so Delta Y has exactly the
    * limits of dy and the pre period still has positive variance.
    args lb sa ia jj n
    clear
    set seed 20260930
    set obs `n'
    gen double id = _n
    gen double z  = 2*rbeta(2,4) - 1
    gen double dy = cond(z >= 0, `ia' + `sa'*z, `lb')
    expand 2
    bysort id: gen byte t = _n - 1
    gen double y = cond(t == 1, dy + rnormal(0, 0.1295), rnormal(0, 0.1295))
    drop dy
end

*-----------------------------------------------------------------------------
* case A: Delta Y^- = -1 exactly (well away from zero), Delta Y^+ = 0.5 + 0.8z
*-----------------------------------------------------------------------------
display as txt _n "== case A: identified set matches the closed form =="
_mkb -1 0.8 0.5 0.5 3000

capture noisily didc_bounds y, runvar(z) time(t) pre(0) post(1) id(id) ///
    h(0.4) b(0.6) c1(0 0.5 1) c2(0 0.5 1)
local rc = _rc
_t didc_bounds_runs 0 `rc'

if `rc' == 0 {
    local mup = e(mu_plus)
    local mum = e(mu_minus)
    local tau = e(tau_didc)
    local bd1 = e(breakdown_c1)
    local bd2 = e(breakdown_c2)
    local ez  = e(identified_set_empty_at_zero)
    display as txt "    mu_plus = " %8.5f `mup' "  mu_minus = " %8.5f `mum' "  tau_didc = " %8.5f `tau'

    * the estimates should be close to the truth (0.5, -1, 1.5)
    local ok = (abs(`mup' - 0.5) < 0.05)
    _t mu_plus_close_to_truth 1 `ok'
    local ok = (abs(`mum' + 1.0) < 0.05)
    _t mu_minus_close_to_truth 1 `ok'
    local ok = (abs(`tau' - 1.5) < 0.05)
    _t tau_close_to_truth 1 `ok'

    matrix LB = e(bounds_lb)
    matrix UB = e(bounds_ub)
    matrix EM = e(bounds_empty)

    local c1v "0 0.5 1"
    local c2v "0 0.5 1"
    local r = 0
    local allok = 1
    foreach a of numlist `c1v' {
        local ++r
        local c = 0
        foreach b of numlist `c2v' {
            local ++c
            local lo = max(`mup' - `a', `tau' - `b')
            local hi = min(`mup' + `a', `tau' + `b')
            if abs(LB[`r',`c'] - `lo') > 1e-10 local allok = 0
            if abs(UB[`r',`c'] - `hi') > 1e-10 local allok = 0
            local expected = (`lo' > `hi')
            if EM[`r',`c'] != `expected' local allok = 0
        }
    }
    _t whole_grid_matches_closed_form 1 `allok'

    local ok = (EM[1,1] == 1)
    _t origin_cell_is_empty 1 `ok'
    local ok = (`ez' == 1)
    _t empty_at_origin_is_reported 1 `ok'
    local ok = (EM[3,3] == 0)
    _t wide_grid_cell_is_not_empty 1 `ok'
    local ok = (abs(`bd1' - `mup') < 1e-10)
    _t breakdown_c1_equals_mu_plus 1 `ok'
    local ok = (abs(`bd2' - `tau') < 1e-10)
    _t breakdown_c2_equals_tau 1 `ok'
}

*-----------------------------------------------------------------------------
* case B: design(rcs) -- the two periods become independent samples
*-----------------------------------------------------------------------------
display as txt _n "== case B: design(rcs) agrees with design(panel) =="
_mkb -1 0.8 0.5 0.5 3000

capture quietly didc y, runvar(z) time(t) pre(0) post(1) id(id) h(0.4) b(0.6)
local rc = _rc
_t panel_runs 0 `rc'
if `rc' == 0 {
    local t_panel  = e(tau_didc)
    local se_panel = e(se_didc)
    local l_panel  = e(lemma1_ok)
    local ok = (`l_panel' == 1)
    _t panel_lemma_check_passes 1 `ok'

    * rebuild with independent samples, then estimate with design(rcs)
    _mkb -1 0.8 0.5 0.5 3000
    replace id = id + 1e6 if t == 1
    capture quietly didc y, runvar(z) time(t) pre(0) post(1) design(rcs) h(0.4) b(0.6)
    local rc2 = _rc
    _t rcs_runs 0 `rc2'
    if `rc2' == 0 {
        local t_rcs  = e(tau_didc)
        local se_rcs = e(se_didc)
        local l_rcs  = e(lemma1_ok)
        display as txt "    panel tau = " %8.5f `t_panel' " (se " %7.5f `se_panel' ")" ///
            "   rcs tau = " %8.5f `t_rcs' " (se " %7.5f `se_rcs' ")"
        local gap = abs(`t_panel' - `t_rcs')
        local tol = 4*(`se_panel' + `se_rcs')
        local ok = (`gap' < `tol')
        _t panel_and_rcs_agree 1 `ok'
        local ok = (`l_rcs' == .)
        _t rcs_lemma_flag_is_missing 1 `ok'
    }
}

display as txt _n "{hline 60}"
display as txt "  passed: " as result "$didc_npass" as txt "    failed: " as result "$didc_nfail"
display as txt "{hline 60}"
if $didc_nfail > 0 {
    display as error "TEST SUITE FAILED"
    exit 9
}
display as result "ALL TESTS PASSED"
