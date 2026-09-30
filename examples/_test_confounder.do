*==============================================================================
* _test_confounder.do -- V2: the whole point of the design.
*
* Under model 1 of the paper there is a TIME-INVARIANT confounder c at the
* cutoff.  Then
*     a plain RD on the post period      estimates  c + tau   (WRONG)
*     a plain RD on the pre period       estimates  c         (the confounder)
*     didc                               estimates  tau       (RIGHT)
* This test asserts all three statements on the same data.
*
*   do "D:/OpenCode/didc/examples/_test_confounder.do"
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

use "D:/OpenCode/didc/data/didc_sim1.dta", clear

local TAU  = TrueTau[1]
local CJS  = TrueConf[1]
display as txt _n "== V2: removing a time-invariant confounder =="
display as txt "    true tau = `TAU'    true confounder c = `CJS'"

*---- didc ---------------------------------------------------------------
capture noisily didc y, runvar(z) time(t) pre(0) post(1) id(id)
local rc = _rc
_t didc_runs 0 `rc'
if `rc' == 0 {
    local t_didc   = e(tau_didc)
    local se_didc  = e(se_didc)
    local mu_plus  = e(mu_plus)
    local mu_minus = e(mu_minus)
    local l1       = e(lemma1_ok)

    display as txt "        didc      tau = " %9.5f `t_didc' "  (se " %7.5f `se_didc' ")"
    display as txt "        mu_plus (Delta Y^+) = " %9.5f `mu_plus'
    display as txt "        mu_minus(Delta Y^-) = " %9.5f `mu_minus'

    local ok = (abs(`t_didc' - `TAU') < 3*`se_didc' + 0.02)
    _t didc_recovers_tau 1 `ok'
    local ok = (`l1' == 1)
    _t lemma1_pass 1 `ok'
}

*---- plain RD on the post period: contaminated by c ---------------------
quietly rdrobust y z if t == 1, c(0)
local rd_post = e(tau_bc)
local se_post = e(se_tau_rb)
display as txt "        RD(y, post only) = " %9.5f `rd_post' "  (se " %7.5f `se_post' ")"
display as txt "        expected c + tau = " %9.5f (`CJS' + `TAU')

local ok = (abs(`rd_post' - (`CJS' + `TAU')) < 3*`se_post' + 0.02)
_t post_only_rd_is_contaminated 1 `ok'

*---- plain RD on the pre period: pure confounder ------------------------
quietly rdrobust y z if t == 0, c(0)
local rd_pre = e(tau_bc)
local se_pre = e(se_tau_rb)
display as txt "        RD(y, pre only)  = " %9.5f `rd_pre' "  (se " %7.5f `se_pre' ")"
display as txt "        expected c       = " %9.5f `CJS'

local ok = (abs(`rd_pre' - `CJS') < 3*`se_pre' + 0.02)
_t pre_only_rd_is_the_confounder 1 `ok'

*---- the gap between them must be the treatment effect ---------------
local gap = `rd_post' - `rd_pre'
local se_gap = sqrt(`se_post'^2 + `se_pre'^2)
display as txt "        RD(post) - RD(pre) = " %9.5f `gap' as txt ///
    "  (expected tau = " %6.3f `TAU' as txt ")"
local ok = (abs(`gap' - `TAU') < 3*`se_gap' + 0.02)
_t difference_of_rds_recovers_tau 1 `ok'

display as txt _n "{hline 60}"
display as txt "  passed: " as result "$didc_npass" as txt "    failed: " as result "$didc_nfail"
display as txt "{hline 60}"
if $didc_nfail > 0 {
    display as error "TEST SUITE FAILED"
    exit 9
}
display as result "ALL TESTS PASSED"
