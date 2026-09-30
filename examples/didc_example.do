*==============================================================================
* didc_example.do -- end-to-end, reproducible demonstration of didc
*
* Run from the package root:
*     do "examples/didc_example.do"
* or from the examples directory:
*     do "didc_example.do"
*
* Covers: the estimate, the Lemma 1 self-check, both validity tests, and the
* partial-identification bounds.  Regenerates its own data if needed.
*==============================================================================

clear all
set more off

*------------------------------------------------------------------------------
* 0. locate the package and the data
*------------------------------------------------------------------------------
capture confirm file "data/didc_sim1.dta"
if _rc {
    capture confirm file "../data/didc_sim1.dta"
    if _rc {
        display as txt "data not found; generating it with examples/didc_simdata.do"
        capture confirm file "examples/didc_simdata.do"
        if _rc == 0 {
            do "examples/didc_simdata.do" 2000 20260930 0.3 0.5 "data"
            local datadir "data"
        }
        else {
            do "didc_simdata.do" 2000 20260930 0.3 0.5 "../data"
            local datadir "../data"
        }
    }
    else {
        local datadir "../data"
    }
}
else {
    local datadir "data"
}

*------------------------------------------------------------------------------
* 1. the setup
*------------------------------------------------------------------------------
* Two periods, t = 0 (pre) and t = 1 (post).  Municipalities (here: units) are
* assigned by z >= 0 to a CONFOUNDER that is present in both periods, and to
* the TREATMENT OF INTEREST that appears only at t = 1.  The confounder adds a
* jump of TrueConf at the cutoff in both periods; the treatment adds TrueTau.

use "`datadir'/didc_sim1.dta", clear
display ""
display as txt "true treatment effect (TrueTau)  = " as result TrueTau[1]
display as txt "true confounding jump (TrueConf) = " as result TrueConf[1]

*------------------------------------------------------------------------------
* 2. what a standard RD would say
*------------------------------------------------------------------------------
* Using only the post-treatment cross-section, the discontinuity is the
* confounder plus the treatment, and the RD is invalid.
display ""
display as txt "{hline 78}"
display as txt "A standard RD on the post-treatment cross-section:"
display as txt "{hline 78}"
rdrobust y z if t == 1, c(0)

*------------------------------------------------------------------------------
* 3. didc
*------------------------------------------------------------------------------
display ""
display as txt "{hline 78}"
display as txt "difference-in-discontinuities:"
display as txt "{hline 78}"
didc y, runvar(z) time(t) pre(0) post(1) id(id) graph(didc_doseplot)

display as txt "e(lemma1_ok) = " as result e(lemma1_ok) as txt ///
    "   (1 = the two parameterisations agree with a common bandwidth)"

display as txt "the graph is stored as didc_doseplot; export it with:"
display as txt "    graph export didc_doseplot.png, replace"

*------------------------------------------------------------------------------
* 4. validity tests (need at least two pre-treatment periods)
*------------------------------------------------------------------------------
use "`datadir'/didc_sim_multi.dta", clear

display ""
display as txt "{hline 78}"
display as txt "is the confounding effect time-invariant?  (stacked-RD Wald test)"
display as txt "{hline 78}"
didc_test y, runvar(z) time(t) test(wald) pre(-1 0)

display ""
display as txt "{hline 78}"
display as txt "is the conditional mean time-invariant?  (bootstrap KS test)"
display as txt "{hline 78}"
didc_test y, runvar(z) time(t) test(ks) pre(-1 0) reps(999)

*------------------------------------------------------------------------------
* 5. partial identification
*------------------------------------------------------------------------------
use "`datadir'/didc_sim1.dta", clear

display ""
didc_bounds y, runvar(z) time(t) pre(0) post(1) id(id) c1(0(0.5)3) c2(0(0.5)3)

display as txt _n "breakdown values: c1 can rise to " as result %6.3f e(breakdown_c1) ///
    as txt " and c2 to " as result %6.3f e(breakdown_c2)
display as txt "before the sign of the effect is no longer identified."
display as txt "the identified set at c1 = c2 = 0 is empty here (" ///
    as result e(identified_set_empty_at_zero) as txt ///
    "), which says the data reject time-invariant confounding."

display as txt _n "didc_example.do: done"
