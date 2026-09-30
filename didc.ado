*! didc 0.1.0  2026-09-30
*! Difference-in-Discontinuities (DiDC) for Stata
*! Picchetti, Pinto & Shinoki (2026), arXiv:2405.18531
*!
*! Point estimation and robust inference are delegated to rdrobust
*! (Calonico, Cattaneo & Titiunik 2014), exactly as in the source paper.
*! What this command adds is the DiDC design layer:
*!   * the differenced-outcome construction (panel) or the common-bandwidth
*!     two-period comparison (repeated cross-sections),
*!   * validation of the structure the design requires,
*!   * an exact check of Lemma 1 (the RD of the differences equals the
*!     difference of the RDs) using a common bandwidth, and
*!   * reporting in DiDC terms: mu_minus = Delta Y^- and mu_plus = Delta Y^+.

program define didc, eclass
    version 17
    * The body runs inside preserve/restore so that a failure cannot leave the
    * user's dataset replaced by the estimation subsample.  Stata has no
    * try/finally, hence the thin wrapper.
    preserve
    capture noisily _didc_body `0'
    local rc = _rc
    capture restore
    if `rc' exit `rc'
end


program define _didc_body, eclass
    version 17

    syntax varname [if] [in] ,                         ///
        RUNvar(varname) TIME(varname)                   ///
        PRE(numlist max=1) POST(numlist max=1)          ///
        [ DESIGN(string) ID(varname) CUToff(real 0)     ///
          P(integer 1) Q(integer 2)                     ///
          KERNEL(string) BWSELECT(string)               ///
          H(numlist min=1 max=2) B(numlist min=1 max=2) ///
          LEVEL(integer 95) GRAPH(string) NOLEMMA1 ENGINE(string) ]

    local depvar `varlist'
    if "`design'" == "" local design "panel"
    local design = lower("`design'")
    if "`engine'" == "" local engine "rdrobust"
    local engine = lower("`engine'")
    if !inlist("`engine'", "rdrobust", "mata") {
        display as error "didc: engine() must be rdrobust (the default) or mata"
        exit 198
    }
    if `q' <= `p' {
        display as error "didc: q() must exceed p() (rdrobust needs q > p)"
        exit 198
    }

    *---- if/in, applied to everything that follows -----------------------
    * built by hand rather than marksample so that observations with a
    * missing outcome are reported by _didc_prep (as dropped units) instead
    * of silently breaking the panel structure check
    tempvar touse
    gen byte `touse' = 1 `if' `in'
    replace `touse' = 0 if missing(`touse')

    preserve
    quietly keep if `touse'

    *---- validate and construct ------------------------------------------
    * optional varname options are passed only when non-empty: id() with an
    * empty value is not accepted by syntax's varname type
    local idopt ""
    if "`id'" != "" local idopt "id(`id')"

    _didc_prep check, depvar(`depvar') runvar(`runvar') time(`time') ///
        pre(`pre') post(`post') design(`design') `idopt' cutoff(`cutoff')

    local n_units = r(n_units)
    local n_below = r(n_below)
    local n_above = r(n_above)

    _didc_prep build, depvar(`depvar') runvar(`runvar') time(`time') ///
        pre(`pre') post(`post') design(`design') `idopt' cutoff(`cutoff')

    *---- estimation -------------------------------------------------------
    local lemma_ok   = .
    local t_delta    = .
    local warn_msg   ""
    local bias_above = .
    local bias_below = .

    if "`design'" == "panel" {

        *---- bandwidths -------------------------------------------------
        * with engine(rdrobust) they come out of the same call that
        * estimates.  With engine(mata) they come from the user, or, when the
        * user does not supply them, from rdrobust used for SELECTION ONLY,
        * so that the built-in engine needs rdrobust just to choose a
        * bandwidth and not to produce a number.
        local eh_l "`h'"
        local eb_l "`b'"
        local sel_only 0
        if "`engine'" == "mata" & ("`h'" == "" | "`b'" == "") {
            capture _didc_lpoly _didc_dy _didc_zc, cutoff(0) p(`p') q(`q') ///
                kernel(`kernel') bwselect(`bwselect') h(`h') b(`b') level(`level')
            if _rc {
                display as error "didc: engine(mata) needs h() and b(), because the"
                display as error "      bandwidth selector could not run (rdrobust missing?)"
                exit 111
            }
            if "`h'" == "" local eh_l "`=r(h_l)' `=r(h_r)'"
            if "`b'" == "" local eb_l "`=r(b_l)' `=r(b_r)'"
            local sel_only 1
        }

        if "`engine'" == "mata" {
            _didc_mata _didc_dy _didc_zc, cutoff(0) p(`p') q(`q') ///
                kernel(`kernel') h(`eh_l') b(`eb_l') level(`level')
        }
        else {
            _didc_lpoly _didc_dy _didc_zc, cutoff(0) p(`p') q(`q')        ///
                kernel(`kernel') bwselect(`bwselect') h(`h') b(`b') level(`level')
        }

        local tau_cl  = r(tau_cl)
        local tau_bc  = r(tau_bc)
        local se_rb   = r(se_rb)
        local ci_l_rb = r(ci_l_rb)
        local ci_r_rb = r(ci_r_rb)
        local mu_plus = r(mu_plus)
        local mu_minus= r(mu_minus)
        local bw_h_l = r(h_l)
        local bw_h_r = r(h_r)
        local bw_b_l = r(b_l)
        local bw_b_r = r(b_r)
        local n_h_l  = r(n_h_l)
        local n_h_r  = r(n_h_r)
        local bias_above = r(bias_above)
        local bias_below = r(bias_below)

        *---- Lemma 1 self-check ------------------------------------------
        * the identity holds at the population level; in finite samples the
        * two parameterisations agree exactly only if a COMMON bandwidth is
        * used, which is why h and b are pinned to the primary run
        if "`nolemma1'" == "" {
            capture which rdrobust
            if _rc {
                display as txt "note: the Lemma 1 cross-check needs rdrobust and was skipped;"
                display as txt "      the built-in engine is validated in examples/_test_engine.do"
            }
            else {
            _didc_lpoly _didc_dy _didc_zc, cutoff(0) p(`p') q(`q') ///
                kernel(`kernel') h(`bw_h_l' `bw_h_r') b(`bw_b_l' `bw_b_r') level(`level')
            local t_dy = r(tau_bc)

            _didc_lpoly _didc_post _didc_zc, cutoff(0) p(`p') q(`q') ///
                kernel(`kernel') h(`bw_h_l' `bw_h_r') b(`bw_b_l' `bw_b_r') level(`level')
            local t_post = r(tau_bc)

            _didc_lpoly _didc_pre _didc_zc, cutoff(0) p(`p') q(`q') ///
                kernel(`kernel') h(`bw_h_l' `bw_h_r') b(`bw_b_l' `bw_b_r') level(`level')
            local t_pre = r(tau_bc)

            local t_delta = `t_post' - `t_pre'
            local gap = abs(`t_dy' - `t_delta')
            if `gap' < 1e-8 {
                local lemma_ok = 1
            }
            else {
                local lemma_ok = 0
                local warn_msg "Lemma 1 check failed: RD of differences = " ///
                    %10.6f `t_dy' " vs difference of RDs = " %10.6f `t_delta'
            }
            }
        }
    }

    else {
        *---- repeated cross-sections: two RDs, one common bandwidth -------
        * Lemma 1 requires a single bandwidth for both periods.  The v1
        * convention is to take it from the post period (where the treatment
        * is active) unless the user pins h()/b() down.
        if "`h'" == "" | "`b'" == "" {
            _didc_lpoly `depvar' _didc_zc if `time' == (`post'), cutoff(0) ///
                p(`p') q(`q') kernel(`kernel') bwselect(`bwselect') ///
                h(`h') b(`b') level(`level')
            if "`h'" == "" {
                local bw_h_l = r(h_l)
                local bw_h_r = r(h_r)
            }
            if "`b'" == "" {
                local bw_b_l = r(b_l)
                local bw_b_r = r(b_r)
            }
        }
        if "`h'" != "" {
            tokenize `h'
            local bw_h_l `1'
            local bw_h_r = cond("`2'" == "", "`1'", "`2'")
        }
        if "`b'" != "" {
            tokenize `b'
            local bw_b_l `1'
            local bw_b_r = cond("`2'" == "", "`1'", "`2'")
        }

        _didc_lpoly `depvar' _didc_zc if `time' == (`post'), cutoff(0) ///
            p(`p') q(`q') kernel(`kernel')                       ///
            h(`bw_h_l' `bw_h_r') b(`bw_b_l' `bw_b_r') level(`level')
        local tau_cl_post  = r(tau_cl)
        local tau_bc_post  = r(tau_bc)
        local se_bc_post   = r(se_rb)
        local mu_plus      = r(mu_plus)

        _didc_lpoly `depvar' _didc_zc if `time' == (`pre'), cutoff(0) ///
            p(`p') q(`q') kernel(`kernel')                       ///
            h(`bw_h_l' `bw_h_r') b(`bw_b_l' `bw_b_r') level(`level')
        local tau_cl_pre   = r(tau_cl)
        local tau_bc_pre   = r(tau_bc)
        local se_bc_pre    = r(se_rb)
        local mu_minus     = r(mu_minus)

        * the two period samples are independent by construction
        local tau_cl  = `tau_cl_post' - `tau_cl_pre'
        local tau_bc  = `tau_bc_post' - `tau_bc_pre'
        local se_rb   = sqrt(`se_bc_post'^2 + `se_bc_pre'^2)
        local zcrit   = invnormal(1 - (100-`level')/200)
        local ci_l_rb = `tau_bc' - `zcrit'*`se_rb'
        local ci_r_rb = `tau_bc' + `zcrit'*`se_rb'

        local t_delta = .
    }

    *---- put the primary results back into e() --------------------------
    if "`design'" == "panel" {
        if "`engine'" == "mata" {
            quietly _didc_mata _didc_dy _didc_zc, cutoff(0) p(`p') q(`q') ///
                kernel(`kernel') h(`bw_h_l' `bw_h_r') b(`bw_b_l' `bw_b_r') level(`level')
        }
        else {
            _didc_lpoly _didc_dy _didc_zc, cutoff(0) p(`p') q(`q') ///
                kernel(`kernel') bwselect(`bwselect') h(`h') b(`b') level(`level')
        }
    }
    else {
        _didc_lpoly `depvar' _didc_zc if `time' == (`post'), cutoff(0) ///
            p(`p') q(`q') kernel(`kernel') h(`bw_h_l' `bw_h_r') ///
            b(`bw_b_l' `bw_b_r') level(`level')
        local n_h_l = e(N_h_l)
        local n_h_r = e(N_h_r)
    }

    *---- augment e() with the DiDC quantities ---------------------------
    ereturn scalar tau_didc    = `tau_cl'
    ereturn scalar tau_didc_bc = `tau_bc'
    ereturn scalar se_didc     = `se_rb'
    ereturn scalar ci_didc_l   = `ci_l_rb'
    ereturn scalar ci_didc_r   = `ci_r_rb'
    ereturn scalar mu_plus     = `mu_plus'
    ereturn scalar mu_minus    = `mu_minus'
    ereturn scalar cutoff      = `cutoff'
    ereturn scalar lemma1_ok   = `lemma_ok'
    ereturn scalar tau_delta_of_rds = `t_delta'
    ereturn scalar n_units     = `n_units'
    ereturn scalar n_below     = `n_below'
    ereturn scalar n_above     = `n_above'
    ereturn scalar bw_h_l      = `bw_h_l'
    ereturn scalar bw_h_r      = `bw_h_r'
    ereturn scalar bw_b_l      = `bw_b_l'
    ereturn scalar bw_b_r      = `bw_b_r'
    ereturn scalar n_h_l       = `n_h_l'
    ereturn scalar n_h_r       = `n_h_r'
    ereturn local  design      "`design'"
    ereturn local  engine      "`engine'"
    ereturn local  kernel      "`kernel'"
    ereturn scalar p_used      = `p'
    ereturn scalar q_used      = `q'
    ereturn scalar level_used  = `level'
    ereturn scalar bias_above  = `bias_above'
    ereturn scalar bias_below  = `bias_below'
    ereturn local  cmd         "didc"
    ereturn local  cmdline     "didc `depvar' `if' `in', runvar(`runvar') time(`time') pre(`pre') post(`post')"

    *---- report ----------------------------------------------------------
    _didc_display

    if "`warn_msg'" != "" {
        display as error "warning: `warn_msg'"
        display as error "         the two parameterisations should agree exactly with a"
        display as error "         common bandwidth; please report this as a bug"
    }

    *---- optional RD plot of the differenced outcome ---------------------
    if "`graph'" != "" {
        if "`design'" == "panel" {
            capture noisily rdplot _didc_dy _didc_zc, c(0) p(`p') ///
                graph_options(name(`graph', replace))
            if _rc {
                display as txt "note: rdplot() could not draw the graph (rc = " _rc ")"
            }
        }
        else {
            display as txt "note: graph() is only available with design(panel)"
        }
    }
end
