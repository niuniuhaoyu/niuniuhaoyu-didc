*==============================================================================
* didc_simdata.do -- simulate data for didc (paper's Monte Carlo models 1-4)
*
* Source: Picchetti, Pinto & Shinoki (2026), arXiv:2405.18531, Appendix D.2
*   n = 1000 (paper), Z ~ 2*Beta(2,4)-1, eps ~ N(0, sigma^2), sigma = 0.1295
*   Y_it = mu_it(Z) + eps_it ,  t = 0 (pre) / 1 (post)
*
* Models
*   1  identical functional forms, TIME-INVARIANT CONFOUNDER c at the cutoff
*   2  identical functional forms, NO confounder
*   3  TIME-VARYING functional forms, time-invariant confounder c
*   4  time-varying functional forms, no confounder
*
* Usage
*   do didc_simdata.do [n] [seed] [tau] [cjs] [path]
*   defaults: n(2000) seed(20260930) tau(0.3) cjs(0.5) path("data")
*
* Output (long panel, two rows per unit):
*   id t z y model TrueTau TrueConf
*   data/didc_sim1.dta ... data/didc_sim4.dta  (+ didc_sim.dta = model 1)
*==============================================================================

program define _didc_sim_one
    * args: modelno n seed tau cjs outfile
    args m n seed tau cjs outfile

    quietly {
        clear
        set seed `seed'
        set obs `n'

        gen double id = _n
        gen double z  = 2*rbeta(2,4) - 1

        * common left-side (z < 0) polynomial, identical in all four models
        gen double beta_z = 0.48 + 1.27*z - 0.5*7.18*z^2 + 0.7*20.21*z^3 ///
                          + 1.1*21.54*z^4 + 1.5*7.33*z^5

        * right-side base polynomial for models 1 and 2
        gen double R1 = 0.52 + 0.84*z - 0.1*3*z^2 - 0.3*7.99*z^3 ///
                      - 0.1*9.01*z^4 + 3.56*z^5

        gen double mu0 = .
        gen double mu1 = .
        replace mu0 = beta_z if z < 0
        replace mu1 = beta_z if z < 0

        if `m' == 1 {
            replace mu0 = R1 + `cjs'         if z >= 0
            replace mu1 = R1 + `cjs' + `tau' if z >= 0
        }
        else if `m' == 2 {
            replace mu0 = R1                 if z >= 0
            replace mu1 = R1 + `tau'         if z >= 0
        }
        else if `m' == 3 {
            replace mu0 = R1                            if z >= 0
            replace mu1 = 0.52 + 0.1*z + `cjs' + `tau'  if z >= 0
        }
        else if `m' == 4 {
            replace mu0 = R1                            if z >= 0
            replace mu1 = 0.52 + 0.1*z + `tau'          if z >= 0
        }
        else {
            display as error "_didc_sim_one: unknown model `m'"
            exit 198
        }

        gen double eps0 = rnormal(0, 0.1295)
        gen double eps1 = rnormal(0, 0.1295)

        * long panel: two rows per unit
        gen byte t = 0
        gen double y = mu0 + eps0
        tempfile pre
        save `pre', replace

        replace t = 1
        replace y = mu1 + eps1
        append using `pre'

        gen byte   model    = `m'
        gen double TrueTau  = `tau'
        gen double TrueConf = cond(`m' == 2 | `m' == 4, 0, `cjs')

        order id t z y model TrueTau TrueConf
        sort id t
        label var z        "running variable (time-invariant)"
        label var y        "outcome"
        label var t        "0 = pre, 1 = post"
        label var model    "1..4 as in Appendix D.2"

        save "`outfile'", replace
    }
    display "wrote `outfile'   N = " _N
end

*------------------------------------------------------------------ main ------
* Positional arguments (Stata's `do` does not accept ", option()"):
*   do didc_simdata.do [n] [seed] [tau] [cjs] [path]
args n seed tau cjs path

if "`n'"    == "" local n    2000
if "`seed'" == "" local seed 20260930
if "`tau'"  == "" local tau  0.3
if "`cjs'"  == "" local cjs  0.5
if "`path'" == "" local path "data"

capture mkdir "`path'"
if _rc {
    * mkdir fails if the directory already exists; that is fine, but a real
    * failure must be reported rather than silently writing nowhere.
    capture dir "`path'"
    if _rc {
        display as error "didc_simdata.do: cannot create or read directory `path'"
        exit 603
    }
}

forvalues m = 1/4 {
    _didc_sim_one `m' `n' `seed' `tau' `cjs' "`path'/didc_sim`m'.dta"
}

* convenience copy: model 1 is the default example dataset
copy "`path'/didc_sim1.dta" "`path'/didc_sim.dta", replace

*------------------------------------------------------------------ multi ------
* Three periods, used by the validity tests, which need at least two
* pre-treatment periods.  t = -1 and t = 0 are pre-treatment and share the
* time-invariant confounder of model 1; the treatment effect is added at
* t = 1.  The functional form above the cutoff is identical in the two pre
* periods, so both validity tests should fail to reject on this dataset.

clear
set seed `seed'
set obs `n'
gen double id = _n
gen double z  = 2*rbeta(2,4) - 1
gen double beta_z = 0.48 + 1.27*z - 0.5*7.18*z^2 + 0.7*20.21*z^3 ///
                  + 1.1*21.54*z^4 + 1.5*7.33*z^5
gen double R1 = 0.52 + 0.84*z - 0.1*3*z^2 - 0.3*7.99*z^3 ///
              - 0.1*9.01*z^4 + 3.56*z^5
gen double mu = .
gen double y  = .
gen byte   t  = .
tempfile m1 m0
replace t  = -1
replace mu = cond(z < 0, beta_z, R1 + `cjs')
replace y  = mu + rnormal(0, 0.1295)
save `m1', replace
replace t  = 0
replace mu = cond(z < 0, beta_z, R1 + `cjs')
replace y  = mu + rnormal(0, 0.1295)
save `m0', replace
replace t  = 1
replace mu = cond(z < 0, beta_z, R1 + `cjs' + `tau')
replace y  = mu + rnormal(0, 0.1295)
append using `m0'
append using `m1'
gen double TrueTau  = `tau'
gen double TrueConf = `cjs'
order id t z y TrueTau TrueConf
sort id t
label var TrueConf "time-invariant confounding jump at the cutoff"
save "`path'/didc_sim_multi.dta", replace

display "didc_simdata.do: done (n=`n', seed=`seed', tau=`tau', cjs=`cjs')"
display "  two-period models 1-4 : didc_sim1..4.dta, didc_sim.dta"
display "  three-period (tests)  : didc_sim_multi.dta"
