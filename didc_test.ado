*! didc_test 0.1.0  2026-09-30
*! Validity tests for the difference-in-discontinuities design.
*! Picchetti, Pinto & Shinoki (2026), arXiv:2405.18531, Section 4.
*!
*!   test(wald)  Section 4.1 -- is the confounding effect time-invariant?
*!               stacked local-linear RDs over the pre-treatment periods,
*!               joint Wald test that the discontinuity is the same in every
*!               pre period.  Rejection means the DiDC design is not credible.
*!
*!   test(ks)    Section 4.2 -- are the conditional mean functions
*!               time-invariant?  Two-sample Kolmogorov-Smirnov statistic with
*!               a bootstrap that imposes the null.  See _didc_ks.ado.
*!
*! Specification actually estimated for test(wald), on |Z - c| <= h, weighted
*! by the kernel:
*!
*!   y = sum_k [ alpha_k + gamma_k D + beta_k (Z-c) + theta_k D (Z-c) ] T_k
*!
*! The paper prints equation (7) without the level term gamma_k D, which would
*! leave only slope changes at the cutoff.  A regression-discontinuity test for
*! a stable discontinuity needs the level jumps, so both are estimated and both
*! are reported: gamma_k (the RD effect in period k) is the headline test and
*! theta_k (the slope change) is reported alongside.  See
*! docs/research-notes.md, section 9, item 8.

program define didc_test, rclass
    version 17
    * A rclass program that returns nothing of its own has its r() -- including
    * that inherited from the subprogram -- cleared on exit.  The body's results
    * must therefore be copied out explicitly before this wrapper returns.
    * Note the asymmetry with didc and didc_bounds: e() survives their wrappers.
    preserve
    capture noisily _didc_test_body `0'
    local rc = _rc
    local vtest "`r(test)'"

    if "`vtest'" == "wald" {
        local vF        = r(F)
        local vdf       = r(df)
        local vp        = r(p)
        local vFt       = r(F_theta)
        local vdft      = r(df_theta)
        local vpt       = r(p_theta)
        local vFj       = r(F_joint)
        local vdfj      = r(df_joint)
        local vpj       = r(p_joint)
        local vbw       = r(bw)
        local vN        = r(N)
        local vkernel   "`r(kernel)'"
        local vperiods  "`r(periods)'"
        tempname mG mT
        capture matrix `mG' = r(gamma)
        capture matrix `mT' = r(theta)

        capture restore
        if `rc' exit `rc'

        return scalar F        = `vF'
        return scalar df       = `vdf'
        return scalar p        = `vp'
        return scalar F_theta  = `vFt'
        return scalar df_theta = `vdft'
        return scalar p_theta  = `vpt'
        return scalar F_joint  = `vFj'
        return scalar df_joint = `vdfj'
        return scalar p_joint  = `vpj'
        return scalar bw       = `vbw'
        return scalar N        = `vN'
        capture return matrix gamma = `mG'
        capture return matrix theta = `mT'
        return local  kernel  "`vkernel'"
        return local  periods "`vperiods'"
        return local  test    "wald"
        exit
    }

    if "`vtest'" == "ks" {
        local vka = r(ks_above)
        local vkb = r(ks_below)
        local vpa = r(p_above)
        local vpb = r(p_below)
        local vna = r(N_above)
        local vnb = r(N_below)
        local vdeg = r(degree)
        local vreps = r(reps)
        local vpn = r(period_new)
        local vpo = r(period_old)

        capture restore
        if `rc' exit `rc'

        return scalar ks_above = `vka'
        return scalar ks_below = `vkb'
        return scalar p_above  = `vpa'
        return scalar p_below  = `vpb'
        return scalar N_above  = `vna'
        return scalar N_below  = `vnb'
        return scalar degree   = `vdeg'
        return scalar reps     = `vreps'
        return scalar period_new = `vpn'
        return scalar period_old = `vpo'
        return local  test     "ks"
        exit
    }

    * the body failed before returning anything
    capture restore
    if `rc' exit `rc'
end


program define _didc_test_body, rclass
    version 17

    syntax varname [if] [in] ,                          ///
        RUNvar(varname) TIME(varname) TEST(string)       ///
        PRE(numlist min=2)                               ///
        [ CUToff(real 0) P(integer 1) KERNEL(string)     ///
          BW(real -1) DEGree(integer 3)                  ///
          REPS(integer 999) SEED(integer 20260930)       ///
          LEVEL(integer 95) NOISily ]

    local depvar `varlist'
    local test = lower("`test'")
    if !inlist("`test'", "wald", "ks") {
        display as error "didc_test: test() must be wald or ks"
        exit 198
    }

    capture which rdrobust
    if _rc {
        display as error "{bf:rdrobust} is required by didc_test but is not installed."
        display as error "    install it with:   ssc install rdrobust, replace"
        exit 111
    }

    *---- the if/in sample ------------------------------------------------
    tempvar touse
    gen byte `touse' = 1 `if' `in'
    replace `touse' = 0 if missing(`touse')

    preserve
    quietly keep if `touse'

    tempvar zc
    gen double `zc' = `runvar' - (`cutoff')

    *---- KS branch -------------------------------------------------------
    if "`test'" == "ks" {
        _didc_ks, depvar(`depvar') runvar(`runvar') time(`time') pre(`pre') ///
            cutoff(`cutoff') degree(`degree') reps(`reps') seed(`seed')
        local ksA = r(ks_above)
        local ksB = r(ks_below)
        local pA  = r(p_above)
        local pB  = r(p_below)
        local nA  = r(N_above)
        local nB  = r(N_below)
        local pn  = r(period_new)
        local po  = r(period_old)

        display ""
        display as txt "Kolmogorov-Smirnov test of time-invariant conditional means"
        display as txt "  periods compared: t = `pn'  vs  t = `po'" ///
            "     basis degree = `degree'     bootstrap reps = `reps'"
        display as txt "{hline 74}"
        display as txt "                     |    KS statistic      p-value      N"
        display as txt "---------------------+-----------------------------------------------"
        display as txt "  above the cutoff   |   " as result %11.6f `ksA' as txt ///
            "   " as result %11.4f `pA' as txt "   " as result %8.0f `nA'
        display as txt "  below the cutoff   |   " as result %11.6f `ksB' as txt ///
            "   " as result %11.4f `pB' as txt "   " as result %8.0f `nB'
        display as txt "{hline 74}"
        display as txt "  H0: the conditional mean of the outcome does not change between"
        display as txt "      the two pre-treatment periods.  A small p-value means the"
        display as txt "      shape of E[Y|Z] moved over time, which undermines the DiDC."

        return scalar ks_above = `ksA'
        return scalar ks_below = `ksB'
        return scalar p_above  = `pA'
        return scalar p_below  = `pB'
        return scalar N_above  = `nA'
        return scalar N_below  = `nB'
        return scalar degree   = `degree'
        return scalar reps     = `reps'
        return scalar period_new = `pn'
        return scalar period_old = `po'
        return local  test "ks"
        exit
    }

    *---- wald branch -----------------------------------------------------
    * common bandwidth: the smallest period-specific CCT bandwidth, unless the
    * user pins it down (Appendix D.1 concludes that any data-driven choice
    * appropriate to that RD works well)
    if `bw' <= 0 {
        local h = 1e300
        foreach v of numlist `pre' {
            capture quietly rdrobust `depvar' `zc' if `time' == `v', c(0) p(`p')
            if _rc == 0 {
                local hv = min(e(h_l), e(h_r))
                if `hv' < `h' local h = `hv'
            }
        }
        if `h' >= 1e300 {
            display as error "didc_test: could not compute a bandwidth; supply bw()"
            exit 2000
        }
    }
    else {
        local h = `bw'
    }

    if "`kernel'" == "" local kernel "triangular"
    local kernel = lower("`kernel'")
    if !inlist("`kernel'", "triangular", "uniform", "epanechnikov") {
        display as error "didc_test: kernel() must be triangular, uniform or epanechnikov"
        exit 198
    }

    * membership in pre() must be built as an explicit indicator: inlist()
    * takes comma-separated arguments, so inlist(t, -1 0) is silently wrong
    tempvar inpre
    gen byte `inpre' = 0
    foreach v of numlist `pre' {
        quietly replace `inpre' = 1 if `time' == `v'
    }
    quietly keep if `inpre' == 1 & abs(`zc') <= `h'

    tempvar w
    if "`kernel'" == "uniform" {
        gen double `w' = 1
    }
    else if "`kernel'" == "triangular" {
        gen double `w' = 1 - abs(`zc')/`h'
    }
    else {
        gen double `w' = 0.75*(1 - (`zc'/`h')^2)
    }

    *---- build the stacked design ---------------------------------------
    local K : word count `pre'
    local i = 0
    local regs ""
    foreach v of numlist `pre' {
        local ++i
        quietly gen double T`i'  = (`time' == `v')
        quietly gen double D`i'  = T`i' * (`zc' >= 0)
        quietly gen double ZL`i' = T`i' * `zc'
        quietly gen double DZ`i' = T`i' * (`zc' >= 0) * `zc'
        local regs "`regs' T`i' D`i' ZL`i' DZ`i'"
    }

    * every period needs observations on both sides of the cutoff
    forvalues i = 1/`K' {
        quietly count if T`i' == 1 & `zc' >= 0
        local nabove = r(N)
        quietly count if T`i' == 1 & `zc' <  0
        local nbelow = r(N)
        if `nabove' < 4 | `nbelow' < 4 {
            display as error "didc_test: period " word("`pre'", `i') ///
                " has too few observations within the bandwidth (`nabove' above, `nbelow' below)"
            exit 2000
        }
    }

    if "`noisily'" != "" {
        regress `depvar' `regs' [aw = `w'], noconstant
    }
    else {
        quietly regress `depvar' `regs' [aw = `w'], noconstant
    }

    *---- the Wald tests --------------------------------------------------
    * Stata's aweights variance treats the kernel weights as if they were
    * inverse-variance weights, i.e. as if Var(u_i) = sigma2 / w_i.  Under the
    * homoskedastic errors of this design that understates sigma2 by the factor
    * sum(w)/n and makes the test reject far too often (measured size about 14%
    * at a nominal 5%).  The kernel enters the ESTIMATOR as a fixed weight, so
    * the correct variance is sigma2 (X'WX)^-1 with
    * sigma2 = sum(w u^2) / (sum(w) - k).  It is built here from X'WX directly.
    quietly predict double _didc_u, residuals
    quietly generate double _didc_wu2 = `w' * _didc_u^2
    quietly summarize _didc_wu2
    local sum_wu2 = r(sum)
    quietly summarize `w'
    local sum_w = r(sum)
    local kpar = `K' * 4
    local dfres = `sum_w' - `kpar'
    if `dfres' <= 0 {
        display as error "didc_test: too few effective observations for `kpar' parameters"
        exit 2000
    }
    local s2 = `sum_wu2' / `dfres'

    matrix accum XWX = `regs' [iw = `w'], noconstant
    matrix V = `s2' * invsym(XWX)
    matrix bvec = e(b)'

    local gtest ""
    local stest ""
    forvalues i = 1/`=`K'-1' {
        local j = `i' + 1
        local gtest "`gtest' (D`i' = D`j')"
        local stest "`stest' (DZ`i' = DZ`j')"
    }

    * R matrices: one row per restriction, with +1 and -1 on the coefficients
    local nr = `K' - 1
    matrix RG = J(`nr', `kpar', 0)
    matrix RS = J(`nr', `kpar', 0)
    forvalues i = 1/`nr' {
        local j = `i' + 1
        matrix RG[`i', `=4*(`i'-1)+2'] =  1
        matrix RG[`i', `=4*(`j'-1)+2'] = -1
        matrix RS[`i', `=4*(`i'-1)+4'] =  1
        matrix RS[`i', `=4*(`j'-1)+4'] = -1
    }

    foreach spec in G S J {
        if "`spec'" == "G" {
            matrix R = RG
        }
        else if "`spec'" == "S" {
            matrix R = RS
        }
        else {
            matrix R = RG \ RS
        }
        matrix Rb = R * bvec
        matrix RVR = R * V * R'
        matrix W = Rb' * invsym(RVR) * Rb
        local F = W[1,1] / rowsof(R)
        local dfn = rowsof(R)
        local p = Ftail(`dfn', `dfres', `F')
        if "`spec'" == "G" {
            local Fg = `F'
            local dg = `dfn'
            local pg = `p'
        }
        else if "`spec'" == "S" {
            local Fs = `F'
            local ds = `dfn'
            local ps = `p'
        }
        else {
            local Fj = `F'
            local dj = `dfn'
            local pj = `p'
        }
    }
    capture drop _didc_u _didc_wu2

    matrix gvec = J(1, `K', .)
    matrix svec = J(1, `K', .)
    forvalues i = 1/`K' {
        matrix gvec[1,`i'] = _b[D`i']
        matrix svec[1,`i'] = _b[DZ`i']
    }
    matrix colnames gvec = `pre'
    matrix colnames svec = `pre'
    matrix rownames gvec = gamma_level_jump
    matrix rownames svec = theta_slope_change

    local hh = `h'
    local lev = `level'

    display ""
    display as txt "Stacked-RD test of a time-invariant discontinuity"
    display as txt "  pre periods: `pre'    bandwidth h = " as result %8.5f `hh' ///
        as txt "    kernel = `kernel'    N = " as result %8.0f _N
    display as txt "{hline 74}"
    display as txt "                       gamma_k (level jump)     theta_k (slope change)"
    display as txt "---------------------+-------------------------------------------------"
    forvalues i = 1/`K' {
        local v : word `i' of `pre'
        display as txt "  t = " as result %-14.0f `v' as txt " |     " ///
            as result %12.6f gvec[1,`i'] as txt "             " ///
            as result %12.6f svec[1,`i']
    }
    display as txt "{hline 74}"
    display as txt "  H0: gamma_1 = ... = gamma_K    F(" as result `dg' as txt "," ///
        as result %4.0f `dfres' as txt ") = " as result %8.4f `Fg' as txt ///
        "    p = " as result %8.4f `pg'
    display as txt "  H0: theta_1 = ... = theta_K    F(" as result `ds' as txt "," ///
        as result %4.0f `dfres' as txt ") = " as result %8.4f `Fs' as txt ///
        "    p = " as result %8.4f `ps'
    display as txt "  H0: both of the above          F(" as result `dj' as txt "," ///
        as result %4.0f `dfres' as txt ") = " as result %8.4f `Fj' as txt ///
        "    p = " as result %8.4f `pj'
    display as txt "{hline 74}"
    display as txt "  Rejection means the discontinuity moved across pre-treatment periods,"
    display as txt "  i.e. the confounding effect is not time-invariant (Assumption 4)."

    return scalar F        = `Fg'
    return scalar df       = `dg'
    return scalar p        = `pg'
    return scalar F_theta  = `Fs'
    return scalar df_theta = `ds'
    return scalar p_theta  = `ps'
    return scalar F_joint  = `Fj'
    return scalar df_joint = `dj'
    return scalar p_joint  = `pj'
    return scalar bw       = `hh'
    return scalar N        = _N
    return matrix gamma = gvec
    return matrix theta = svec
    return local  kernel "`kernel'"
    return local  periods "`pre'"
    return local  test "wald"
end
