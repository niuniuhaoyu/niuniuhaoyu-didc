*! _didc_display 0.1.0  2026-09-30
*! Internal: print the didc result table.  Reads everything from e().

program define _didc_display
    version 17

    local design   "`e(design)'"
    local lev      = e(level_used)
    local kern     "`e(kernel)'"
    local pord     = e(p_used)
    local qord     = e(q_used)
    local tau      = e(tau_didc)
    local taubc    = e(tau_didc_bc)
    local se       = e(se_didc)
    local cil      = e(ci_didc_l)
    local cir      = e(ci_didc_r)
    local mup      = e(mu_plus)
    local mum      = e(mu_minus)
    local c        = e(cutoff)
    local lemma_ok = e(lemma1_ok)
    local tdelta   = e(tau_delta_of_rds)
    local n_units  = e(n_units)
    local n_l      = e(n_below)
    local n_r      = e(n_above)
    local nhl      = e(n_h_l)
    local nhr      = e(n_h_r)
    local bwhl     = e(bw_h_l)
    local bwhr     = e(bw_h_r)
    local bbl      = e(bw_b_l)
    local bbr      = e(bw_b_r)

    local z = `taubc' / `se'
    local pval = 2*normal(-abs(`z'))
    local fmt "%9.4f"

    if "`design'" == "panel" {
        local what "RD of the differences  [design(panel)]"
    }
    else {
        local what "difference of the RDs   [design(rcs)]"
    }

    display ""
    display as txt "Difference-in-discontinuities estimation"
    display as txt "Estimator: `what'"
    display as txt "Cutoff c = " as result %8.4g `c' as txt ///
        "    Kernel = " as result "`kern'" as txt ///
        "    Order est. (p) = " as result "`pord'" as txt ///
        "    Order bias (q) = " as result "`qord'" as txt ///
        "    Engine = " as result "`e(engine)'"
    if "`design'" == "panel" {
        display as txt "Units  = " as result %8.0g `n_units' as txt ///
            "    below c = " as result %8.0g `n_l' as txt ///
            "    above c = " as result %8.0g `n_r'
    }
    else {
        display as txt "Observations below c = " as result %8.0g `n_l' as txt ///
            "    above c = " as result %8.0g `n_r'
    }
    display as txt "BW est. (h): below " as result %8.4g `bwhl' as txt ///
        "  above " as result %8.4g `bwhr' as txt ///
        "     BW bias (b): below " as result %8.4g `bbl' as txt ///
        "  above " as result %8.4g `bbr'
    display as txt "Eff. obs within h: below " as result %8.0g `nhl' as txt ///
        "  above " as result %8.0g `nhr'

    display as txt "{hline 78}"
    display as txt "                          |  Estimate   Robust Std.Err.      z     P>|z|"
    display as txt "                          |                       " ///
        "[" as result "`lev'%" as txt " Conf. Interval]"
    display as txt "--------------------------+-----------------------------------------------"
    display as txt "  tau_didc (conventional) | " as result `fmt' `tau'
    display as txt "  tau_didc (bias-corrected)| " as result `fmt' `taubc' as txt ///
        "      " as result `fmt' `se' as txt "   " ///
        as result %8.2f `z' as txt "   " as result %7.3f `pval'
    display as txt "                          |                       " ///
        "[" as result `fmt' `cil' as txt ", " as result `fmt' `cir' as txt "]"
    display as txt "{hline 78}"
    display as txt "  mu_minus = Delta Y^- (below c) = " as result `fmt' `mum'
    display as txt "  mu_plus  = Delta Y^+ (above c) = " as result `fmt' `mup'
    display as txt "  note: tau_didc = mu_plus - mu_minus by construction"

    if "`lemma_ok'" == "." {
        display as txt "  Lemma 1 check: not available (design(rcs) has no Delta Y)"
    }
    else if `lemma_ok' == 1 {
        display as txt "  Lemma 1 check: " as result "PASS" as txt ///
            " (RD of differences = difference of RDs = " as result `fmt' `tdelta' as txt ")"
    }
    else {
        display as error "  Lemma 1 check: FAIL"
    }
    display as txt "{hline 78}"
end
