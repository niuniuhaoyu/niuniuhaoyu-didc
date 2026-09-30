*==============================================================================
* _test_vs_rdrobust.do -- V3, R leg
*
* Exports the simulated data so that examples/reference/run_rdrobust.R can
* compare didc's numbers with the R package rdrobust.  R is not required by
* the package and is not installed on every machine, so this file only does
* the export; the comparison that needs no R lives in
* examples/_test_vs_wls.do, which checks the same conventional local-linear
* estimate against a hand-rolled weighted least squares fit.
*
*   do "D:/OpenCode/didc/examples/_test_vs_rdrobust.do"
*   then:  cd reference && Rscript run_rdrobust.R
*==============================================================================

clear all
set more off
adopath + "D:/OpenCode/didc"

capture confirm file "D:/OpenCode/didc/data/didc_sim1.dta"
if _rc {
    display as error "_test_vs_rdrobust.do: data/didc_sim1.dta not found;"
    display as error "    run examples/didc_simdata.do first"
    exit 601
}

use "D:/OpenCode/didc/data/didc_sim1.dta", clear
keep id t z y
export delimited id t z y using "D:/OpenCode/didc/examples/reference/didc_sim1.csv", replace

display as txt "wrote examples/reference/didc_sim1.csv"
display as txt "now run:  cd examples/reference && Rscript run_rdrobust.R"

*------------------------------------------------------------------------------
* Same comparison without R (the in-repository check).
*------------------------------------------------------------------------------
display as txt _n "running the R-free equivalent, examples/_test_vs_wls.do"
do "D:/OpenCode/didc/examples/_test_vs_wls.do"
