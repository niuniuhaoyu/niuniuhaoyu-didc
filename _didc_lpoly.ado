*! _didc_lpoly 0.1.0  2026-09-30
*! Internal: run one local-polynomial RD on (outcome, running variable) and
*! return the pieces didc needs, including the left and right intercepts of
*! the CONVENTIONAL local polynomial (which are Delta Y^- and Delta Y^+ in
*! the notation of Picchetti, Pinto & Shinoki 2026, and which the
*! partial-identification bounds are built from).
*!
*! All estimation is delegated to rdrobust (Calonico, Cattaneo & Titiunik
*! 2014).  The error message is rewritten so that failures are attributable
*! to didc rather than to rdrobust.

program define _didc_lpoly, rclass
    version 17

    syntax varlist(min=2 max=2 numeric) [if] [in], ///
        [ CUToff(real 0) P(integer 1) Q(integer 2) KERNEL(string) ///
          BWSELECT(string) H(numlist min=1 max=2) B(numlist min=1 max=2) ///
          LEVEL(integer 95) ]

    tokenize `varlist'
    local dv `1'
    local rv `2'

    local opts "c(`cutoff') p(`p') q(`q') level(`level')"
    if "`kernel'"   != "" local opts "`opts' kernel(`kernel')"
    if "`bwselect'" != "" local opts "`opts' bwselect(`bwselect')"
    if "`h'"        != "" local opts "`opts' h(`h')"
    if "`b'"        != "" local opts "`opts' b(`b')"

    capture quietly rdrobust `dv' `rv' `if' `in', `opts'
    if _rc {
        local rc = _rc
        display as error "didc: local polynomial estimation failed (rc = `rc')"
        if `rc' == 2000 | `rc' == 111 | `rc' == 148 {
            display as error "     usually too few observations in the bandwidth"
            display as error "     around c(`cutoff'); try a larger h() or a smaller p()"
        }
        exit `rc'
    }

    *---- conventional (not bias-corrected) results -----------------------
    return scalar tau_cl  = e(tau_cl)
    return scalar se_cl   = e(se_tau_cl)
    return scalar ci_l_cl = e(ci_l_cl)
    return scalar ci_r_cl = e(ci_r_cl)

    *---- bias-corrected, robust inference --------------------------------
    return scalar tau_bc  = e(tau_bc)
    return scalar se_rb   = e(se_tau_rb)
    return scalar ci_l_rb = e(ci_l_rb)
    return scalar ci_r_rb = e(ci_r_rb)
    return scalar pv_rb   = e(pv_rb)

    *---- bandwidths and effective sample sizes ---------------------------
    return scalar h_l   = e(h_l)
    return scalar h_r   = e(h_r)
    return scalar b_l   = e(b_l)
    return scalar b_r   = e(b_r)
    return scalar n_h_l = e(N_h_l)
    return scalar n_h_r = e(N_h_r)
    return scalar n_l   = e(N_l)
    return scalar n_r   = e(N_r)

    *---- left/right intercepts of the conventional fit -------------------
    * these are mu_minus = Delta Y^- and mu_plus = Delta Y^+ of the paper
    tempname BPr BPl
    matrix `BPr' = e(beta_Y_p_r)
    matrix `BPl' = e(beta_Y_p_l)
    return scalar mu_plus  = `BPr'[1,1]
    return scalar mu_minus = `BPl'[1,1]

    *---- the estimated bias on each side, for cross-engine checks -------
    return scalar bias_above = e(bias_r)
    return scalar bias_below = e(bias_l)

    return scalar p = e(p)
    return scalar q = e(q)
    return local  kernel   "`e(kernel)'"
    return local  bwselect "`e(bwselect)'"
    return scalar level = e(level)
end
