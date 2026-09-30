*! _didc_ks 0.1.0  2026-09-30
*! Internal: two-sample Kolmogorov-Smirnov test of time-invariance of the
*! conditional mean functions (Section 4.2 of Picchetti, Pinto & Shinoki 2026).
*!
*! Statistic, above the cutoff (below is symmetric):
*!
*!   KS_+ = sup_{z in S+(Z)} (1/sqrt(n+)) | sum_{i: z0 <= Z_i <= z}
*!                                        ( p(Z_i)'ghat_0 - p(Z_i)'ghat_-1 ) |
*!
*! with p(.) a polynomial basis of degree deg and ghat_t the least squares
*! estimate of the conditional mean in period t.  The statistic is therefore
*! the largest absolute cumulative sum of d_i = p(Z_i)'(ghat_0 - ghat_-1),
*! with the observations of that side ordered by Z.
*!
*! The p-value comes from a bootstrap that IMPOSES the null of a common
*! conditional mean: the pooled two-period data on each side estimate the
*! common mean, residuals are resampled with replacement, and the two period
*! means are re-estimated in every bootstrap sample.
*!
*! Implementation notes
*!   * The printed KS_- in the paper sums over Z_i >= z0, contradicting its own
*!     definition of S^-(Z).  We use Z_i < z0, symmetric with KS_+.  See
*!     docs/research-notes.md, section 9, item 2.
*!   * In Mata, X[k, ] with k a 0/1 vector selects rows BY POSITION.  Logical
*!     subsetting must go through select().  This file uses select() throughout.

version 17

mata:
mata clear

real matrix _didc_ks_basis(real colvector zz, real scalar deg)
{
    real matrix X
    real scalar j
    X = J(rows(zz), deg + 1, .)
    for (j = 1; j <= deg + 1; j++) X[, j] = zz:^(j - 1)
    return(X)
}

real scalar _didc_ks_stat(real colvector zz, real colvector dd)
{
    real colvector ord, ds, cs
    ord = order(zz, 1)
    ds  = dd[ord]
    cs  = runningsum(ds)
    return(max(abs(cs)))
}

void _didc_ks_run(string scalar yvar, string scalar zvar, string scalar tvar,
                  real scalar deg, real scalar p1, real scalar p2,
                  real scalar reps, real scalar seed)
{
    real colvector y, zc, t, sel, mA, mB
    real colvector zA, yA, tA, zB, yB, tB
    real scalar nA, nB, b, ksA, ksB, ksAstar, ksBstar, pA, pB
    real matrix XA, XB, XA0, XAm1, XB0, XBm1
    real matrix iA0, iAm1, iB0, iBm1
    real colvector gA0, gAm1, gB0, gBm1, gAb, gBb
    real colvector uA, uB, d, ustar, yAstar, yBstar, h0, hm1, idx
    real colvector yA0, yAm1, yB0, yBm1, yA0s, yAm1s, yB0s, yBm1s

    y  = st_data(., yvar)
    zc = st_data(., zvar)
    t  = st_data(., tvar)

    sel = (t :== p1) :| (t :== p2)
    y = select(y, sel)
    zc = select(zc, sel)
    t = select(t, sel)

    mA = (zc :>= 0)
    mB = (zc :<  0)
    zA = select(zc, mA); yA = select(y, mA); tA = select(t, mA)
    zB = select(zc, mB); yB = select(y, mB); tB = select(t, mB)

    nA = rows(zA)
    nB = rows(zB)

    if (nA < 6 * (deg + 1) | nB < 6 * (deg + 1)) {
        errprintf("_didc_ks: too few observations on one side of the cutoff\n")
        _error(2000)
    }

    XA = _didc_ks_basis(zA, deg)
    XB = _didc_ks_basis(zB, deg)

    XA0  = select(XA, (tA :== p1)); yA0  = select(yA, (tA :== p1))
    XAm1 = select(XA, (tA :== p2)); yAm1 = select(yA, (tA :== p2))
    XB0  = select(XB, (tB :== p1)); yB0  = select(yB, (tB :== p1))
    XBm1 = select(XB, (tB :== p2)); yBm1 = select(yB, (tB :== p2))

    if (rows(yA0) < deg + 2 | rows(yAm1) < deg + 2 | rows(yB0) < deg + 2 | rows(yBm1) < deg + 2) {
        errprintf("_didc_ks: one period has too few observations on one side\n")
        _error(2000)
    }

    iA0  = invsym(XA0'XA0)
    iAm1 = invsym(XAm1'XAm1)
    iB0  = invsym(XB0'XB0)
    iBm1 = invsym(XBm1'XBm1)

    gA0  = iA0  * (XA0'yA0)
    gAm1 = iAm1 * (XAm1'yAm1)
    gB0  = iB0  * (XB0'yB0)
    gBm1 = iBm1 * (XBm1'yBm1)

    d   = XA * (gA0 - gAm1)
    ksA = _didc_ks_stat(zA, d) / sqrt(nA)
    d   = XB * (gB0 - gBm1)
    ksB = _didc_ks_stat(zB, d) / sqrt(nB)

    /* pooled fits impose the null of a common conditional mean */
    gAb = invsym(XA'XA) * (XA'yA)
    gBb = invsym(XB'XB) * (XB'yB)
    uA = yA - XA * gAb
    uB = yB - XB * gBb

    rseed(seed)
    pA = 0
    pB = 0
    for (b = 1; b <= reps; b++) {
        ustar  = uA[ceil(nA * runiform(nA, 1))]
        yAstar = XA * gAb + ustar

        ustar  = uB[ceil(nB * runiform(nB, 1))]
        yBstar = XB * gBb + ustar

        yA0s  = select(yAstar, (tA :== p1))
        yAm1s = select(yAstar, (tA :== p2))
        yB0s  = select(yBstar, (tB :== p1))
        yBm1s = select(yBstar, (tB :== p2))

        h0  = iA0  * (XA0'yA0s)
        hm1 = iAm1 * (XAm1'yAm1s)
        d   = XA * (h0 - hm1)
        ksAstar = _didc_ks_stat(zA, d) / sqrt(nA)
        if (ksAstar >= ksA) pA = pA + 1

        h0  = iB0  * (XB0'yB0s)
        hm1 = iBm1 * (XBm1'yBm1s)
        d   = XB * (h0 - hm1)
        ksBstar = _didc_ks_stat(zB, d) / sqrt(nB)
        if (ksBstar >= ksB) pB = pB + 1
    }

    st_numscalar("_didc_ks_a",  ksA)
    st_numscalar("_didc_ks_b",  ksB)
    st_numscalar("_didc_ks_pa", (pA + 1) / (reps + 1))
    st_numscalar("_didc_ks_pb", (pB + 1) / (reps + 1))
    st_numscalar("_didc_ks_na", nA)
    st_numscalar("_didc_ks_nb", nB)
}
end


program define _didc_ks, rclass
    version 17

    syntax , depvar(varname) runvar(varname) time(varname) ///
             pre(numlist min=2) [ CUToff(real 0) DEGree(integer 3) ///
             REPS(integer 999) SEED(integer 20260930) ]

    * use the two most recent pre-treatment periods: the test compares the
    * period just before the treatment with the one before that
    local np : word count `pre'
    local last  = -1e300
    local prev  = -1e300
    foreach v of numlist `pre' {
        if `v' > `last' {
            local prev `last'
            local last `v'
        }
        else if `v' > `prev' {
            local prev `v'
        }
    }
    if `prev' == -1e300 {
        display as error "didc_test: pre() must supply at least two distinct periods"
        exit 198
    }

    quietly count if `time' == `last'
    if r(N) == 0 {
        display as error "didc_test: no observations with `time' == `last'"
        exit 2000
    }
    quietly count if `time' == `prev'
    if r(N) == 0 {
        display as error "didc_test: no observations with `time' == `prev'"
        exit 2000
    }

    tempvar zc touse
    gen double `zc' = `runvar' - (`cutoff')
    gen byte `touse' = 1
    replace `touse' = 0 if missing(`depvar') | missing(`runvar') | missing(`time')

    preserve
    quietly keep if `touse'
    capture noisily mata: _didc_ks_run("`depvar'", "`zc'", "`time'", `degree', `last', `prev', `reps', `seed')
    local rc = _rc
    capture restore
    if `rc' exit `rc'

    local ksA = scalar(_didc_ks_a)
    local ksB = scalar(_didc_ks_b)
    local pA  = scalar(_didc_ks_pa)
    local pB  = scalar(_didc_ks_pb)
    local nA  = scalar(_didc_ks_na)
    local nB  = scalar(_didc_ks_nb)

    return scalar ks_above  = `ksA'
    return scalar ks_below  = `ksB'
    return scalar p_above   = `pA'
    return scalar p_below   = `pB'
    return scalar N_above   = `nA'
    return scalar N_below   = `nB'
    return scalar degree    = `degree'
    return scalar reps      = `reps'
    return scalar period_new = `last'
    return scalar period_old = `prev'
end
