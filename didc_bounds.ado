*! didc_bounds 0.1.0  2026-09-30
*! Partial identification and sensitivity analysis for the
*! difference-in-discontinuities design.
*! Picchetti, Pinto & Shinoki (2026), arXiv:2405.18531, Section 5.1.
*!
*! Assumption 4 (time-invariant confounding) is replaced by bounded variation:
*!
*!   A6  | E[Y_1(1,0) - Y_0(1,0) | Z=z] | <= c1          (a time trend in the
*!                                                        confounded outcome)
*!       | [confounding effect at t=1] - [confounding effect at t=0] | <= c2
*!                                                        (drift in the
*!                                                        confounding effect)
*!
*! Lemma 6 then gives the identified set
*!
*!   tau_c in [ max{ mu_plus - c1 , tau_didc - c2 } ,
*!              min{ mu_plus + c1 , tau_didc + c2 } ]
*!
*! with mu_plus = Delta Y^+ and tau_didc = Delta Y^+ - Delta Y^-.
*!
*! The CONVENTIONAL (not bias-corrected) local polynomial estimates are used,
*! because those are the ones in equation (4) of the paper that Lemma 6 plugs
*! into max/min. See docs/research-notes.md, section 11.2.
*!
*! An EMPTY identified set (lower bound above upper bound) is a real output,
*! not an error: it means the data are incompatible with that pair of
*! sensitivity parameters (equivalently, that Assumption 4 is rejected, which
*! is exactly what the Wald test of didc_test checks). The algebraic condition
*! is Delta Y^- > 0 at c1 = c2 = 0. See research-notes.md, section 11.4.

program define didc_bounds, eclass
    version 17
    * see the note in didc.ado about this wrapper
    preserve
    capture noisily _didc_bounds_body `0'
    local rc = _rc
    capture restore
    if `rc' exit `rc'
end


program define _didc_bounds_body, eclass
    version 17

    syntax varname [if] [in] ,                            ///
        RUNvar(varname) TIME(varname)                      ///
        PRE(numlist max=1) POST(numlist max=1)             ///
        [ DESIGN(string) ID(varname) CUToff(real 0)        ///
          C1(numlist) C2(numlist)                          ///
          P(integer 1) Q(integer 2)                        ///
          KERNEL(string) BWSELECT(string)                  ///
          H(numlist min=1 max=2) B(numlist min=1 max=2)    ///
          LEVEL(integer 95) MATRIX(name) ]

    local depvar `varlist'
    if "`c1'" == "" local c1 "0(0.5)5"
    if "`c2'" == "" local c2 "0(0.5)5"

    *---- the DiDC estimate, reusing didc itself --------------------------
    local idopt ""
    if "`id'" != "" local idopt "id(`id')"
    local hop ""
    if "`h'" != "" local hop "h(`h')"
    local bop ""
    if "`b'" != "" local bop "b(`b')"
    local kop ""
    if "`kernel'" != "" local kop "kernel(`kernel')"
    local bwop ""
    if "`bwselect'" != "" local bwop "bwselect(`bwselect')"

    * nolemma1: the self check costs three extra regressions and adds nothing
    * to the bounds, which inherit whatever the primary estimate is
    quietly didc `depvar' `if' `in', runvar(`runvar') time(`time') ///
        pre(`pre') post(`post') design(`design') `idopt' ///
        cutoff(`cutoff') p(`p') q(`q') `kop' `bwop' `hop' `bop' ///
        level(`level') nolemma1

    local mu_plus  = e(mu_plus)
    local mu_minus = e(mu_minus)
    local tau      = e(tau_didc)

    local n1 : word count `c1'
    local n2 : word count `c2'

    matrix LB = J(`n1', `n2', .)
    matrix UB = J(`n1', `n2', .)
    matrix EM = J(`n1', `n2', .)
    local rlab ""
    local clab ""

    local i = 0
    foreach a of numlist `c1' {
        local ++i
        local rlab "`rlab' `a'"
        local j = 0
        foreach b of numlist `c2' {
            local ++j
            local lo = max(`mu_plus' - `a', `tau' - `b')
            local hi = min(`mu_plus' + `a', `tau' + `b')
            matrix LB[`i',`j'] = `lo'
            matrix UB[`i',`j'] = `hi'
            matrix EM[`i',`j'] = (`lo' > `hi')
        }
    }
    local j = 0
    foreach b of numlist `c2' {
        local ++j
        local clab "`clab' `b'"
    }
    matrix rownames LB = `rlab'
    matrix colnames LB = `clab'
    matrix rownames UB = `rlab'
    matrix colnames UB = `clab'
    matrix rownames EM = `rlab'
    matrix colnames EM = `clab'

    * breakdown values: the largest sensitivity parameters at which the
    * conclusion "the effect is positive" still holds, holding the other
    * parameter at zero.  LB > 0 requires mu_plus - c1 > 0 AND tau - c2 > 0.
    local bd_c1 = `mu_plus'
    local bd_c2 = `tau'

    * At c1 = c2 = 0 the identified set is the single point tau_c and is EMPTY
    * whenever Delta Y^- is not exactly zero: max{mu_plus,tau} > min{mu_plus,tau}
    * iff mu_plus != tau iff mu_minus != 0.  The sign of mu_minus is irrelevant.
    local empty_now = ( max(`mu_plus', `tau') > min(`mu_plus', `tau') )

    display ""
    display as txt "Partial identification under bounded variation  (Section 5.1, Lemma 6)"
    display as txt "  mu_plus = Delta Y^+ = " as result %10.6f `mu_plus' ///
        as txt "     mu_minus = Delta Y^- = " as result %10.6f `mu_minus'
    display as txt "  tau_didc = mu_plus - mu_minus = " as result %10.6f `tau'
    display as txt "  c1 bounds a time trend in the confounded outcome;"
    display as txt "  c2 bounds drift in the confounding effect itself."
    display as txt "{hline 82}"

    local j = 0
    display as txt "   identified set [LB, UB] as a function of c1 (rows) and c2"
    display as txt "   (columns).  Each cell shows LB , UB ;  an 'empty' flag means"
    display as txt "   LB > UB, i.e. the data are incompatible with that (c1, c2)."
    display as txt "{hline 82}"
    display as txt "   c1 \ c2   " _continue
    foreach b of numlist `c2' {
        display as txt %9.2f `b' _continue
    }
    display ""
    local i = 0
    foreach a of numlist `c1' {
        local ++i
        display as txt "  " %7.2f `a' "   " _continue
        local j = 0
        foreach b of numlist `c2' {
            local ++j
            if EM[`i',`j'] == 1 {
                display as error %9.2f LB[`i',`j'] as txt "," _continue
            }
            else {
                display as result %9.2f LB[`i',`j'] as txt "," _continue
            }
        }
        display ""
    }
    display as txt "{hline 82}"

    * a compact second table: the lower bound only, which is what decides
    * whether a qualitative conclusion survives
    display as txt "  lower bound LB(c1, c2) = max{ mu_plus - c1 , tau_didc - c2 }"
    display as txt "   c1 \ c2   " _continue
    foreach b of numlist `c2' {
        display as txt %9.2f `b' _continue
    }
    display ""
    local i = 0
    foreach a of numlist `c1' {
        local ++i
        display as txt "  " %7.2f `a' "   " _continue
        local j = 0
        foreach b of numlist `c2' {
            local ++j
            display as result %9.2f LB[`i',`j'] _continue
        }
        display ""
    }
    display as txt "{hline 82}"
    display as txt "  breakdown values for the conclusion 'tau_c > 0':"
    display as txt "    c1 can rise to " as result %10.6f `bd_c1' as txt " (holding c2 = 0)"
    display as txt "    c2 can rise to " as result %10.6f `bd_c2' as txt " (holding c1 = 0)"
    if `empty_now' == 1 {
        display as error "  note: at c1 = c2 = 0 the identified set is EMPTY"
        display as error "        (mu_minus = Delta Y^- = " %8.4f `mu_minus' " is not exactly zero),"
        display as error "        which is direct evidence that Assumption 4 fails in these data."
    }

    if "`matrix'" != "" {
        matrix `matrix' = LB
    }

    ereturn matrix bounds_lb = LB
    ereturn matrix bounds_ub = UB
    ereturn matrix bounds_empty = EM
    ereturn scalar mu_plus  = `mu_plus'
    ereturn scalar mu_minus = `mu_minus'
    ereturn scalar tau_didc = `tau'
    ereturn scalar breakdown_c1 = `bd_c1'
    ereturn scalar breakdown_c2 = `bd_c2'
    ereturn scalar identified_set_empty_at_zero = `empty_now'
    ereturn local  cmd "didc_bounds"
end
