*! trop: Triply Robust Panel Estimators
*! Version 0.1.1 May 19, 2026
*! Author: Clarke Damian, Justin Waddy
*! dclarke@fen.uchile.cl, j.waddy@exeter.ac.uk

/*
Versions
0.1.0 May 05, 2025: Pre-release
*/

capture program drop trop
program define trop, eclass
    version 15.0

    #delimit ;
    syntax varlist(min=4 max=4 numeric) [if] [in],
        [
        lambda_unit(string)
        lambda_time(string)
        lambda_nn(string)
        treated_periods(integer 10)
        vce(string)
        reps(integer 200)
        seed(integer 0)
        cv(string)
        unit_grid(numlist)
        time_grid(numlist)
        nn_grid(numlist)
        n_treated_units(integer 1)
        kfold(integer 0)
        level(integer 95)
        returnweights
        generate(string)
        verbose
        ]
        ;
    #delimit cr

    *--------------------------------------------------------------------------*
    * (0) Parse inputs
    *--------------------------------------------------------------------------*
    tokenize `varlist'
    local outcome  `1'
    local unitvar  `2'
    local timevar  `3'
    local treatvar `4'

    tempvar touse
    mark `touse' `if' `in'

    *--------------------------------------------------------------------------*
    * (1) Work only on estimation sample
    *--------------------------------------------------------------------------*
    preserve
        quietly keep if `touse'

        * If unit ID is string, convert to numeric group ID
        capture confirm numeric variable `unitvar'
        if _rc {
            tempvar _unit_id
            egen `_unit_id' = group(`unitvar')
            local unitvar `_unit_id'
        }

        *------------------------------------------------------------------*
        * (2) Basic panel and variable checks
        *------------------------------------------------------------------*
        quietly xtset `unitvar' `timevar'

        if `"`r(balanced)'"' != "strongly balanced" {
            di as error "Panel is unbalanced. trop currently requires a strongly balanced panel."
            exit 451
        }

        quietly count if missing(`outcome')
        if r(N) != 0 {
            di as error "Missing values found in dependent variable. A balanced panel without missing observations is required."
            exit 416
        }

        quietly count if missing(`treatvar')
        if r(N) != 0 {
            di as error "Missing values found in treatment variable. A balanced panel without missing observations is required."
            exit 416
        }

        quietly count if !inlist(`treatvar', 0, 1)
        if r(N) != 0 {
            di as error "Treatment variable takes values distinct from 0 and 1."
            exit 450
        }

        quietly summarize `treatvar', meanonly
        if r(min) == 0 & r(max) == 0 {
            di as error "All units are controls."
            exit 459
        }
        if r(min) == 1 & r(max) == 1 {
            di as error "All units are treated."
            exit 459
        }

        if `treated_periods' <= 0 {
            di as error "treated_periods() must be a positive integer."
            exit 198
        }

        *------------------------------------------------------------------*
        * (3) Robust panel dimensions
        *------------------------------------------------------------------*
        tempvar unit_tag time_tag
        quietly egen `unit_tag' = tag(`unitvar')
        quietly egen `time_tag' = tag(`timevar')

        quietly count if `unit_tag'
        local N = r(N)

        quietly count if `time_tag'
        local T = r(N)

        drop `unit_tag' `time_tag'

        quietly count
        local NT = r(N)

        if `NT' != `N' * `T' {
            di as error "Internal panel dimension error: observations are not equal to N*T."
            di as error "Detected N = `N', T = `T', observations = `NT'."
            exit 459
        }

        if `treated_periods' >= `T' {
            di as error "treated_periods() must be smaller than the total number of time periods."
            di as error "Detected T = `T', treated_periods = `treated_periods'."
            exit 198
        }

        if "`verbose'" != "" {
            di as txt "Detected N units       = `N'"
            di as txt "Detected T periods     = `T'"
            di as txt "Detected observations  = `NT'"
            di as txt "Treated post periods   = `treated_periods'"
        }

        *------------------------------------------------------------------*
        * (4) Tuning parameters
        *------------------------------------------------------------------*
        local lambda_unit_set = 0
        if "`lambda_unit'" != "" {
            capture confirm number `lambda_unit'
            if _rc {
                di as error "lambda_unit() must be a nonnegative real."
                exit 198
            }
            if `lambda_unit' < 0 {
                di as error "lambda_unit() must be nonnegative."
                exit 198
            }
            local lambda_unit_set = 1
        }

        local lambda_time_set = 0
        if "`lambda_time'" != "" {
            capture confirm number `lambda_time'
            if _rc {
                di as error "lambda_time() must be a nonnegative real."
                exit 198
            }
            if `lambda_time' < 0 {
                di as error "lambda_time() must be nonnegative."
                exit 198
            }
            local lambda_time_set = 1
        }

        * lambda_nn is parsed for now but not yet used by this minimal WLS implementation
        local lambda_nn_set = 0
        local lambda_nn_inf = 0
        if "`lambda_nn'" == "" {
            * unspecified -> later CV placeholder
        }
        else if inlist(lower("`lambda_nn'"), "inf", "infinity", ".") {
            local lambda_nn_inf = 1
            local lambda_nn_set = 1
        }
        else {
            capture confirm number `lambda_nn'
            if _rc {
                di as error "lambda_nn() must be a nonnegative real, or one of {inf, infinity, .}."
                exit 198
            }
            if `lambda_nn' < 0 {
                di as error "lambda_nn() must be nonnegative."
                exit 198
            }
            local lambda_nn_set = 1
        }

        if `lambda_unit_set' local lu = `lambda_unit'
        else                 local lu = 0

        if `lambda_time_set' local lt = `lambda_time'
        else                 local lt = 0

        *------------------------------------------------------------------*
        * (5) Identify treated units
        *------------------------------------------------------------------*
        tempvar ever_treated unit_row unit_tag2
        quietly bysort `unitvar': egen `ever_treated' = max(`treatvar')
        quietly egen `unit_row' = group(`unitvar')
        quietly egen `unit_tag2' = tag(`unitvar')

        quietly levelsof `unit_row' if `ever_treated' == 1 & `unit_tag2' == 1, local(treated_rows)
        local n_treated : word count `treated_rows'

        if `n_treated' == 0 {
            di as error "No treated units found."
            exit 459
        }

        if "`verbose'" != "" {
            di as txt "Detected treated units = `n_treated'"
            di as txt "Treated rows           = `treated_rows'"
        }

        *------------------------------------------------------------------*
        * (6) Wide matrices
        *------------------------------------------------------------------*
        sort `unitvar' `timevar'

        mata: Y = rowshape(st_data(., "`outcome'"), `N')
        mata: W = rowshape(st_data(., "`treatvar'"), `N')

        if "`verbose'" != "" {
            mata: printf("rows(Y), cols(Y) = %g, %g\n", rows(Y), cols(Y))
            mata: printf("rows(W), cols(W) = %g, %g\n", rows(W), cols(W))
        }

        *------------------------------------------------------------------*
        * (7) Estimate
        *------------------------------------------------------------------*
        mata: treated_units = strtoreal(tokens(st_local("treated_rows")))'

        mata: delta = .; delta_unit = .; delta_time = .
        mata: trop_cell_weights(Y, treated_units, `lu', `lt', `treated_periods', ///
                                 delta, delta_unit, delta_time)

        mata: tau_hat = .; mu_hat = .; alpha_hat = .; beta_hat = .
        mata: trop_fit_wls(Y, W, delta, tau_hat, mu_hat, alpha_hat, beta_hat)

        mata: st_numscalar("r(tau_hat)", tau_hat)
        scalar tau_hat_scalar = r(tau_hat)

        * Optional generated fitted residual components can be added later.
        if "`generate'" != "" {
            di as error "generate() is not implemented yet."
            exit 198
        }

        *------------------------------------------------------------------*
        * (8) Return results
        *------------------------------------------------------------------*
        ereturn clear
        ereturn scalar tau = tau_hat_scalar
        ereturn scalar N = `N'
        ereturn scalar T = `T'
        ereturn scalar N_treated = `n_treated'
        ereturn scalar treated_periods = `treated_periods'
        ereturn scalar lambda_unit = `lu'
        ereturn scalar lambda_time = `lt'

        ereturn local cmd "trop"
        ereturn local outcome "`outcome'"
        ereturn local unitvar "`unitvar'"
        ereturn local timevar "`timevar'"
        ereturn local treatvar "`treatvar'"

        di as txt _newline "Triply Robust Panel Estimator"
        di as txt "--------------------------------"
        di as txt "ATT estimate = " as result %9.6f e(tau)

    restore
end


*------------------------------------------------------------------------------*
* Mata functions
*------------------------------------------------------------------------------*

capture mata: mata drop trop_cell_weights()
capture mata: mata drop trop_fit_wls()

mata:

void trop_cell_weights(
    real matrix Y,
    real colvector treated_units,
    real scalar lambda_unit,
    real scalar lambda_time,
    real scalar treated_periods,
    real matrix delta,
    real colvector delta_unit,
    real rowvector delta_time
)
{
    real scalar N, T, T0, center
    real rowvector dist_time, avg_treated
    real colvector dist_unit, A, B
    real matrix mask, sqdiff

    N  = rows(Y)
    T  = cols(Y)
    T0 = T - treated_periods

    if (treated_periods <= 0 | treated_periods >= T) {
        errprintf("treated_periods must be positive and smaller than T.\n")
        exit(198)
    }

    if (rows(treated_units) == 0) {
        errprintf("No treated units supplied to trop_cell_weights().\n")
        exit(459)
    }

    center    = T - treated_periods / 2
    dist_time = abs(((1..T) :- 1) :- center)

    mask = J(N, T, 1)
    mask[., (T0 + 1)..T] = J(N, treated_periods, 0)

    avg_treated = mean(Y[treated_units, .])
    sqdiff    = (J(N, 1, 1) * avg_treated - Y) :^ 2
    A         = rowsum(sqdiff :* mask)
    B         = rowsum(mask)
    dist_unit = sqrt(A :/ B)

    delta_unit = exp(-lambda_unit :* dist_unit)
    delta_time = exp(-lambda_time :* dist_time)

    delta = delta_unit * delta_time
}


void trop_fit_wls(
    real matrix Y,
    real matrix W,
    real matrix delta,
    real scalar tau,
    real scalar mu,
    real colvector alpha,
    real rowvector beta
)
{
    real scalar N, T, NT, p
    real colvector y_vec, w_vec, b
    real matrix X, D_unit, D_time, XtWX, XtWy

    N  = rows(Y)
    T  = cols(Y)
    NT = N * T

    y_vec = vec(Y')
    w_vec = vec(delta')

    D_unit = I(N)[., 2..N] # J(T, 1, 1)
    D_time = J(N, 1, 1) # I(T)[., 2..T]

    p = 1 + (N - 1) + (T - 1) + 1
    X = J(NT, 1, 1), D_unit, D_time, vec(W')

    XtWX = quadcross(X, w_vec, X)
    XtWy = quadcross(X, w_vec, y_vec)

    b = invsym(XtWX) * XtWy

    mu    = b[1]
    alpha = 0 \ b[2..N]
    beta  = (0, b[(N + 1)..(N + T - 1)]')
    tau   = b[p]
}

end
