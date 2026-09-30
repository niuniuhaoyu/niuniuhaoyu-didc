*==============================================================================
* _test_all.do -- run the whole didc test suite in order
*
*   do "D:/OpenCode/didc/examples/_test_all.do" [mcreps]
*
* Each script exits with a nonzero code if an assertion fails, so as a batch
* job this file stops at the first failure and the log shows where.
*==============================================================================

clear all
set more off
adopath + "D:/OpenCode/didc"

local mcreps 500
args a1
if "`a1'" != "" local mcreps `a1'

local here : pwd
capture cd "D:/OpenCode/didc/examples"
if _rc {
    display as error "_test_all.do: cannot find D:/OpenCode/didc/examples"
    exit 601
}

display as txt _n "{hline 78}"
display as txt "  didc test suite"
display as txt "{hline 78}"

foreach f in _test_prep _test_closedform _test_confounder _test_vs_wls ///
             _test_bounds _test_engine _test_validity_size {
    display as txt _n ">>> `f'.do"
    capture noisily do "`f'.do"
    if _rc {
        display as error _n "FAILED: `f'.do (rc = " _rc ")"
        exit _rc
    }
}

display as txt _n ">>> _test_mc_coverage.do `mcreps'"
capture noisily do "_test_mc_coverage.do" `mcreps'
if _rc {
    display as error _n "FAILED: _test_mc_coverage.do (rc = " _rc ")"
    exit _rc
}

display as txt _n "{hline 78}"
display as result "  ENTIRE SUITE PASSED"
display as txt "{hline 78}"
