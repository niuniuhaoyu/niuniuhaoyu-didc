*! _didc_prep 0.1.0  2026-09-30
*! Internal helpers for didc: input validation and construction of the
*! differenced outcome.  Not intended to be called directly by users.
*!
*! Subcommands
*!   _didc_prep check, depvar() runvar() time() pre() post() [design() id() cutoff()]
*!   _didc_prep build, depvar() runvar() time() pre() post() [design() id() cutoff()]
*!
*! Design notes
*!   * panel: the same units are observed in both periods.  Delta Y is formed
*!     within unit and a single RD is run on (Delta Y, Z) -- the left-hand side
*!     of Lemma 1 of Picchetti, Pinto & Shinoki (2026).
*!   * rcs: repeated cross-sections; the two periods are independent samples.
*!     Delta Y cannot be formed within unit, so one RD per period is run with a
*!     common bandwidth and the two are differenced -- the right-hand side of
*!     Lemma 1.
*!
*!   The runvar must be time-invariant (an assumption of the source paper).
*!   check enforces this; build relies on it.
*!
*! Stata gotchas found while developing this file (kept as a warning to future
*! maintainers):
*!   * local x : lower "..." is NOT accepted by Stata 19 (r101);
*!     use local x = lower("...") instead.
*!   * egen nvals() is not available.
*!   * egen ... = tag(...) sets r(N) to 0, not to the number of tagged
*!     observations; always follow it with an explicit count.

program define _didc_prep
    version 17
    gettoken sub 0 : 0, parse(" ,")
    local sub = lower("`sub'")
    if "`sub'" == "check" {
        _didc_check `0'
    }
    else if "`sub'" == "build" {
        _didc_build `0'
    }
    else {
        display as error "_didc_prep: unknown subcommand `sub'"
        exit 198
    }
end


*==============================================================================
* check -- validate the input dataset; abort with an actionable message
*==============================================================================
program define _didc_check, rclass
    version 17

    syntax , depvar(varname) runvar(varname) time(varname) ///
             pre(numlist max=1) post(numlist max=1)        ///
             [ design(string) id(varname) cutoff(real 0) ]

    *---- design and dependency ------------------------------------------
    if "`design'" == "" local design "panel"
    local design = lower("`design'")
    if !inlist("`design'", "panel", "rcs") {
        display as error "didc: design() must be {bf:panel} or {bf:rcs}"
        exit 198
    }
    if "`design'" == "panel" & "`id'" == "" {
        display as error "didc: design(panel) requires id();"
        display as error "      use design(rcs) for repeated cross-sections"
        exit 198
    }
    if "`design'" == "rcs" & "`id'" != "" {
        display as error "didc: id() is not allowed with design(rcs)"
        exit 198
    }

    capture which rdrobust
    if _rc {
        display as error "{bf:rdrobust} is required by didc but is not installed."
        display as error "    install it with:   ssc install rdrobust, replace"
        exit 111
    }

    *---- periods ---------------------------------------------------------
    if (`pre') == (`post') {
        display as error "didc: pre() and post() must be different"
        exit 198
    }

    quietly count if missing(`time')
    if r(N) > 0 {
        display as error "didc: time variable {bf:`time'} has " r(N) " missing value(s)"
        exit 198
    }
    quietly count if `time' != (`pre') & `time' != (`post')
    if r(N) > 0 {
        display as error "didc: {bf:`time'} takes " r(N) " observation(s) outside pre()/post()"
        display as error "     didc supports exactly two periods"
        exit 198
    }

    quietly count if `time' == (`pre')
    local n_pre = r(N)
    quietly count if `time' == (`post')
    local n_post = r(N)
    if `n_pre' == 0 {
        display as error "didc: no observations in the pre period (`time' == `pre')"
        exit 2000
    }
    if `n_post' == 0 {
        display as error "didc: no observations in the post period (`time' == `post')"
        exit 2000
    }

    quietly count if missing(`depvar')
    local n_miss = r(N)
    if `n_miss' > 0 {
        display as txt "note: " `n_miss' " observation(s) have missing {bf:`depvar'}; those units are dropped."
    }

    *---- panel structure --------------------------------------------------
    local n_units = .
    local junk ""
    if "`design'" == "panel" {
        tempvar nper dupt tag1
        local junk "`junk' `nper' `dupt' `tag1'"
        quietly bysort `id': gen long `nper' = _N
        quietly bysort `id' (`time'): gen byte `dupt' = (`time'[1] == `time'[_N])
        quietly egen byte `tag1' = tag(`id') if (`nper' != 2) | `dupt'
        quietly count if `tag1'
        local n_bad = r(N)
        if `n_bad' > 0 {
            display as error "didc: " `n_bad' " unit(s) in {bf:`id'} lack exactly one observation"
            display as error "     in each of the pre and post periods"
            exit 459
        }

        tempvar zmin zmax tag2
        local junk "`junk' `zmin' `zmax' `tag2'"
        quietly bysort `id': egen double `zmin' = min(`runvar')
        quietly bysort `id': egen double `zmax' = max(`runvar')
        quietly egen byte `tag2' = tag(`id') if `zmin' != `zmax'
        quietly count if `tag2'
        local n_tvar = r(N)
        if `n_tvar' > 0 {
            display as error "didc: runvar {bf:`runvar'} varies over time within " `n_tvar' " unit(s)"
            display as error "     the DiDC design requires a time-invariant running variable"
            exit 459
        }

        quietly count
        local n_units = r(N)/2
    }

    *---- support around the cutoff ---------------------------------------
    if "`design'" == "panel" {
        tempvar tag3 tag4
        local junk "`junk' `tag3' `tag4'"
        quietly egen byte `tag3' = tag(`id') if `runvar' <  (`cutoff')
        quietly count if `tag3' == 1
        local n_below = r(N)
        quietly egen byte `tag4' = tag(`id') if `runvar' >= (`cutoff')
        quietly count if `tag4' == 1
        local n_above = r(N)
    }
    else {
        quietly count if `runvar' <  (`cutoff')
        local n_below = r(N)
        quietly count if `runvar' >= (`cutoff')
        local n_above = r(N)
    }

    capture drop `junk'

    if `n_below' == 0 {
        display as error "didc: no observations below cutoff()" 
        exit 2000
    }
    if `n_above' == 0 {
        display as error "didc: no observations above cutoff()"
        exit 2000
    }
    if `n_below' < 10 | `n_above' < 10 {
        display as txt "warning: only " `n_below' " below and " `n_above' " above the cutoff" ///
            " (" cond("`design'" == "panel", "units", "observations") ");" ///
            " estimates may be unreliable."
    }

    return scalar n_pre   = `n_pre'
    return scalar n_post  = `n_post'
    return scalar n_below = `n_below'
    return scalar n_above = `n_above'
    return scalar n_units = `n_units'
    return scalar n_missing_outcome = `n_miss'
    return local  design  "`design'"
end


*==============================================================================
* build -- form the estimation sample in memory
*   panel : one row per unit with _didc_dy = y(post) - y(pre) and
*           _didc_zc = Z - cutoff
*   rcs   : input unchanged plus _didc_zc (the caller runs one RD per period)
* Callers must wrap this in preserve/restore.
*==============================================================================
program define _didc_build, rclass
    version 17

    syntax , depvar(varname) runvar(varname) time(varname) ///
             pre(numlist max=1) post(numlist max=1)        ///
             [ design(string) id(varname) cutoff(real 0) ]

    if "`design'" == "" local design "panel"
    local design = lower("`design'")

    capture drop _didc_dy
    capture drop _didc_zc

    * the running variable is time-invariant, so centering is unambiguous
    gen double _didc_zc = `runvar' - (`cutoff')
    label var _didc_zc "running variable, centered at cutoff()"

    if "`design'" == "rcs" {
        return scalar n_units = .
        return local  design  "rcs"
        exit
    }

    *---- panel: Delta Y within unit, then one row per unit ----------------
    keep if `time' == (`pre') | `time' == (`post')

    tempvar ispost tagm
    gen byte `ispost' = (`time' == (`post'))

    * within unit, sorted by `ispost', the last record is post and the first
    * is pre.  Keep both levels as well as their difference: the difference
    * feeds the estimator, and the two levels feed the Lemma 1 self-check.
    quietly bysort `id' (`ispost'): gen double _didc_pre  = `depvar'[1]
    quietly bysort `id' (`ispost'): gen double _didc_post = `depvar'[_N]
    gen double _didc_dy = _didc_post - _didc_pre

    quietly egen byte `tagm' = tag(`id') if missing(_didc_dy)
    quietly count if `tagm'
    local n_drop = r(N)
    drop if missing(_didc_dy)

    quietly by `id': keep if _n == 1
    keep `id' _didc_zc _didc_pre _didc_post _didc_dy
    quietly count
    local n_units = r(N)

    label var _didc_dy   "y(post) - y(pre)"
    label var _didc_pre  "y in the pre period"
    label var _didc_post "y in the post period"

    return scalar n_units = `n_units'
    return scalar n_dropped_missing_dy = `n_drop'
    return local  design "panel"
end
