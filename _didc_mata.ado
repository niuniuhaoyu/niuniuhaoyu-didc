*! _didc_mata 0.2.0  2026-09-30
*! Internal: the self-contained local polynomial engine for didc.
*!
*! This file computes, without calling rdrobust, everything didc needs from a
*! local polynomial regression of an outcome on a running variable at a cutoff:
*!
*!   1. the conventional left and right intercepts (mu_minus, mu_plus) and the
*!      conventional RD effect,
*!   2. the bias-corrected RD effect, using a q-order fit at bandwidth b,
*!   3. a conventional and a robust standard error.
*!
*! Design
*!   Every quantity is a LINEAR functional of the outcome with known weights.
*!   The point estimates come from weighted least squares; the bias correction
*!   subtracts a known multiple of the q-order derivative estimate; and the
*!   variance is built from the weights themselves, which makes it the EXACT
*!   finite-sample variance of the estimator for a fixed design rather than an
*!   asymptotic approximation:
*!
*!       Var(mu_s)   = s2_s * sum_i w_si^2
*!       Var(tau)    = sum_s Var(mu_s)
*!
*!   Note the difference from rdrobust, which reports CCT's asymptotic variance
*!   and by default estimates sigma^2 by nearest neighbours.  The two agree as
*!   the sample grows but need not agree in finite samples.  Point estimates
*!   and bias corrections, by contrast, reproduce rdrobust exactly; see
*!   examples/_test_engine.do.
*!
*! Side conventions, matching rdrobust
*!   "above" is z >= cutoff, "below" is z < cutoff.  Inside a side the running
*!   variable is oriented so that u = side*(z - cutoff)/bandwidth >= 0, which is
*!   why the below-side derivative carries a factor (-1)^(p+1).

version 17

mata:
mata clear

/* kernel weights, truncated so that they are exactly zero outside |u| <= 1 */
real colvector _didc_K(real colvector u, string scalar kern)
{
    real colvector a, w
    a = abs(u)
    if (kern == "triangular")        w = (1 :- a) :* (a :<= 1)
    else if (kern == "uniform")      w = (a :<= 1)
    else                             w = 0.75 :* (1 :- u:^2) :* (a :<= 1)
    return(w)
}

real matrix _didc_poly(real colvector u, real scalar p)
{
    real matrix X
    real scalar j
    X = J(rows(u), p + 1, 1)
    for (j = 2; j <= p + 1; j++) X[, j] = u:^(j - 1)
    return(X)
}

/*  One side.  side = +1 above the cutoff, -1 below.
 *  sgn maps the u-scale derivative of order p+1 onto the z scale.
 *  Returns, in a row vector:
 *     mu, var_conventional, wss_robust, s2, bias, mu_bc, n
 */
real rowvector _didc_side(real colvector y, real colvector z, real scalar h,
                          real scalar b, real scalar p, real scalar q,
                          real scalar side, real scalar sgn, string scalar kern)
{
    real scalar p1
    p1 = p + 1

    real colvector m, yy, u, ub, k, kb
    m  = (z :* side) :>= 0
    yy = select(y, m)
    u  = (select(z, m) :* side) / h
    ub = (select(z, m) :* side) / b
    k  = _didc_K(u, kern)
    kb = _didc_K(ub, kern)

    real matrix X, Xq, A, AQ, AA, AQQ
    X   = _didc_poly(u, p)
    Xq  = _didc_poly(ub, q)
    A   = quadcross(X, k, X)
    AQ  = quadcross(Xq, kb, Xq)
    AA  = invsym(A)
    AQQ = invsym(AQ)

    real colvector bb, bq
    bb = AA  * (X'  * (k  :* yy))
    bq = AQQ * (Xq' * (kb :* yy))

    real colvector res
    real scalar s2
    res = yy - X * bb
    s2  = sum(k :* res:^2) / sum(k)

    real rowvector e0, eP1
    e0  = J(1, p + 1, 0); e0[1] = 1
    eP1 = J(1, q + 1, 0); eP1[p1 + 1] = 1

    /* conventional weights and their variance factor */
    real rowvector wc
    wc = e0 * AA * (X' * diag(k))
    real scalar varc
    varc = s2 * sum(wc:^2)

    /* CCT bias constant B_s = p! e_0' Gamma_s^-1 vartheta_{s,p+1} */
    real scalar Bconst
    Bconst = factorial(p) * (e0 * AA * (X' * (k :* u:^p1)))

    /* the z-scale (p+1)-th derivative and the resulting bias */
    real scalar deriv, bias
    deriv = sgn * bq[p1 + 1] / b^p1
    bias  = h^p1 * Bconst * deriv

    /* robust weights: the bias-corrected side estimate is linear in y */
    real rowvector wr
    wr = wc - h^p1 * Bconst * sgn * (eP1 * AQQ * (Xq' * diag(kb)))

    return(bb[1], varc, s2 * sum(wr:^2), s2, bias, bb[1] - bias, rows(yy))
}
end

program define _didc_mata, rclass
    version 17

    syntax varlist(min=2 max=2 numeric) [if] [in], ///
        [ CUToff(real 0) P(integer 1) Q(integer 2) KERNEL(string) ///
          H(numlist min=1 max=2) B(numlist min=1 max=2) LEVEL(integer 95) ]

    tokenize `varlist'
    local dv `1'
    local rv `2'

    if "`kernel'" == "" local kernel "triangular"
    local kernel = lower("`kernel'")
    if !inlist("`kernel'", "triangular", "uniform", "epanechnikov") {
        display as error "_didc_mata: kernel() must be triangular, uniform or epanechnikov"
        display as error "      other kernels are available through engine(rdrobust)"
        exit 198
    }
    if "`h'" == "" | "`b'" == "" {
        display as error "_didc_mata: h() and b() are required by the built-in engine"
        exit 198
    }
    if `q' <= `p' {
        display as error "_didc_mata: q() must exceed p()"
        exit 198
    }
    if `p' != 1 {
        display as error "_didc_mata: the built-in engine reproduces rdrobust exactly for"
        display as error "      p(1), which is the default.  For higher-order polynomials"
        display as error "      use engine(rdrobust): the bias correction for p >= 2 is not"
        display as error "      yet bit-for-bit identical (measured gap of about 6e-04 in"
        display as error "      the bias-corrected estimate at p=2, q=3)."
        exit 198
    }

    tokenize `h'
    local hl `1'
    local hr `2'
    if "`hr'" == "" local hr `hl'
    tokenize `b'
    local bl `1'
    local br `2'
    if "`br'" == "" local br `bl'

    * centre the running variable; the two sides are chosen inside Mata
    tempvar zc
    gen double `zc' = `rv' - (`cutoff')

    quietly count if `zc' <  0 & abs(`zc') <= `hl'
    local n_l = r(N)
    quietly count if `zc' >= 0 & abs(`zc') <= `hr'
    local n_r = r(N)

    mata: _didc_mata_run("`dv'", "`zc'", `p', `q', `hl', `hr', `bl', `br', "`kernel'")

    local tau_cl  = scalar(_didc_tau_cl)
    local tau_bc  = scalar(_didc_tau_bc)
    local se_cl   = scalar(_didc_se_cl)
    local se_rb   = scalar(_didc_se_rb)
    local mu_p    = scalar(_didc_mu_plus)
    local mu_m    = scalar(_didc_mu_minus)
    local bias_r  = scalar(_didc_bias_above)
    local bias_l  = scalar(_didc_bias_below)
    local n_r     = scalar(_didc_n_above)
    local n_l     = scalar(_didc_n_below)
    local zcrit   = invnormal(1 - (100-`level')/200)

    quietly count if `zc' <  0 & abs(`zc') <= `hl'
    local nhl = r(N)
    quietly count if `zc' >= 0 & abs(`zc') <= `hr'
    local nhr = r(N)

    return scalar tau_cl  = `tau_cl'
    return scalar tau_bc  = `tau_bc'
    return scalar se_cl   = `se_cl'
    return scalar se_rb   = `se_rb'
    return scalar ci_l_cl = `tau_cl' - `zcrit'*`se_cl'
    return scalar ci_r_cl = `tau_cl' + `zcrit'*`se_cl'
    return scalar ci_l_rb = `tau_bc' - `zcrit'*`se_rb'
    return scalar ci_r_rb = `tau_bc' + `zcrit'*`se_rb'
    return scalar mu_plus  = `mu_p'
    return scalar mu_minus = `mu_m'
    return scalar bias_above = `bias_r'
    return scalar bias_below = `bias_l'
    return scalar h_l = `hl'
    return scalar h_r = `hr'
    return scalar b_l = `bl'
    return scalar b_r = `br'
    return scalar n_l = `n_l'
    return scalar n_r = `n_r'
    return scalar n_h_l = `nhl'
    return scalar n_h_r = `nhr'
    return scalar pv_rb = 2*normal(-abs(`tau_bc'/`se_rb'))
    return scalar level = `level'
    return local  kernel "`kernel'"
    return local  engine "mata"
end

mata:
void _didc_mata_run(string scalar yv, string scalar zv, real scalar p, real scalar q,
                    real scalar hl, real scalar hr, real scalar bl, real scalar br,
                    string scalar kern)
{
    real colvector y, z
    y = st_data(., yv)
    z = st_data(., zv)

    real rowvector R, L
    R = _didc_side(y, z, hr, br, p, q,  1,  1,            kern)
    L = _didc_side(y, z, hl, bl, p, q, -1,  (-1)^(p+1),   kern)

    st_numscalar("_didc_mu_plus",   R[1])
    st_numscalar("_didc_mu_minus",  L[1])
    st_numscalar("_didc_tau_cl",    R[1] - L[1])
    st_numscalar("_didc_tau_bc",    R[6] - L[6])
    st_numscalar("_didc_se_cl",     sqrt(R[2] + L[2]))
    st_numscalar("_didc_se_rb",     sqrt(R[3] + L[3]))
    st_numscalar("_didc_bias_above", R[5])
    st_numscalar("_didc_bias_below", L[5])
    st_numscalar("_didc_s2_above",   R[4])
    st_numscalar("_didc_s2_below",   L[4])
    st_numscalar("_didc_n_above",    R[7])
    st_numscalar("_didc_n_below",    L[7])
}
end
