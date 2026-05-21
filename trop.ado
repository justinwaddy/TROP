*! trop: Triply Robust Panel Estimators
*! Version 0.x.x May 5, 2025 
*! Author: Clarke Damian, Justin Waddy
*! dclarke@fen.uchile.cl, j.waddy@exeter.ac.uk

/*
Versions
0.1.0 May 05, 2025: Pre-release
*/

cap program drop trop
program trop, eclass
version 15.0

#delimit ;
syntax varlist(min=4 max=4) [if] [in],
    [
    lambda_unit(string)
    lambda_time(string)
    lambda_nn(string)
    vce(string)
    reps(integer 200)
    seed(integer 0)
    cv(string)
    unit_grid(numlist)
    time_grid(numlist)
    nn_grid(numlist)
    kfold(integer 0)
    level(integer 95)
    returnweights
    generate(string)
    verbose
    ]
    ;
#delimit cr

*Above we'll follow the same core inputs: outcome variable, unit ID, time variable, and treatment variable.
*Can discuss later whether they explicitly want the lambda inputs to always be chosen by the user

*------------------------------------------------------------------------------*
* (0) Error checks
*------------------------------------------------------------------------------*
*Temporary dataset called touse
tempvar touse
mark `touse' `if' `in'

*Declare boolean, check if numeric, create numberic id if not. 
local stringvar=0 
cap count if `2'==0
if _rc!=0 {
    local stringvar=1
    local groupvar `2'
    tempvar ID
    egen `ID' = group(`2')
    local varlist `1' `ID' `3' `4'
    tokenize `varlist'
}
else {
    tokenize `varlist'
}

if (length("`if'")+length("`in'")>0) {
    preserve
    qui keep if `touse'
}

qui xtset `2' `3'
if `"`r(balanced)'"'!="strongly balanced" {
    dis as error "Panel is unbalanced."
    exit 451
}
tempvar unit_tag time_tag
quietly egen `unit_tag' = tag(`2')
quietly egen `time_tag' = tag(`3')

quietly count if `unit_tag'
local N = r(N)

quietly count if `time_tag'
local T = r(N)

drop `unit_tag' `time_tag'

qui count if `1'==.
if r(N)!=0 {
    dis as error "Missing values found in dependent variable.  A balanced panel without missing observations is required."
    exit 416
}
qui count if `4'==.
if r(N)!=0 {
    dis as error "Missing values found in treatment variable.  A balanced panel without missing observations is required."
    exit 416
}
qui count if `4'!=0&`4'!=1
if r(N)!=0 {
    dis as error "Treatment variable takes values distinct from 0 and 1."
    exit 450
}
qui sum `4'
if (r(min)==0 & r(max)==0)==1 {
    di as error "All units are controls."
    exit 459
}
qui sum `4'
if (r(min)==1 & r(max)==1)==1 {
    di as error "All units are treated."
    exit 459
}

*Restore dataset
if (length("`if'")+length("`in'")>0) {
    restore
}

*--- Tuning parameters ---*
local lambda_unit_set = 0
if "`lambda_unit'" != "" {
    cap confirm number `lambda_unit'
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
    cap confirm number `lambda_time'
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

* lambda_nn: empty (CV chooses), "inf"/"infinity"/"." (drop L), else nonneg real
local lambda_nn_set = 0
local lambda_nn_inf = 0
if "`lambda_nn'" == "" {
    * unspecified -> resolved by CV later
}
else if inlist(lower("`lambda_nn'"), "inf", "infinity", ".") {
    local lambda_nn_inf = 1
    local lambda_nn_set = 1
}
else {
    cap confirm number `lambda_nn'
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



*------------------------------------------------------------------------------*
* (1) Set-up 
*------------------------------------------------------------------------------*
*Prepare dataset before treatment. See example prep: https://github.com/ostasovskyi/TROP-Estimator/blob/main/notebooks/tutorial.ipynb?short_path=b4fd5c7
if (length("`if'")+length("`in'")>0) {
    preserve
    qui keep if `touse'
}

*Wide panels
mata: Y = rowshape(st_data(., "`1'"), `N')
mata: W = rowshape(st_data(., "`4'"), `N')

tempvar ever_treated unit_row unit_tag

qui bys `2': egen `ever_treated' = max(`4')
qui egen `unit_row' = group(`2')
qui egen `unit_tag' = tag(`2')

qui levelsof `unit_row' if `ever_treated' == 1 & `unit_tag' == 1, local(treated_rows)
local n_treated : word count `treated_rows'

drop `ever_treated' `unit_row' `unit_tag'

if (length("`if'")+length("`in'")>0) {
    restore
}

*------------------------------------------------------------------------------*
* (2) Calculate ATT
*------------------------------------------------------------------------------*
*First subsection will be cross-validation code to tune parameters for lambda_unit, lambda_time and lambda_nn if they are not provided.

*CV Placeholder

*Core function will be in this section
if `lambda_unit_set' local lu = `lambda_unit'
else                 local lu = 0

if `lambda_time_set' local lt = `lambda_time'
else                 local lt = 0

mata: delta = .; delta_unit = .; delta_time = .
mata: trop_cell_weights(Y, W, strtoreal(tokens(st_local("treated_rows")))', ///
                        `lu', `lt', delta, delta_unit, delta_time)
mata: tau_hat = .; mu_hat = .; alpha_hat = .; beta_hat = .; L_hat = .; iters = .
if `lambda_nn_inf' {
    mata: trop_fit_wls(Y, W, delta, tau_hat, mu_hat, alpha_hat, beta_hat)
}
else {
    mata: trop_fit_nuclear(Y, W, delta, `lambda_nn', 1e-10, 5000, tau_hat, mu_hat, alpha_hat, beta_hat, L_hat, iters)
}
mata: st_local("tau_hat", strofreal(tau_hat))
*--------------------------------------------------------------------------*
* (3) Standard error: bootstrap
*--------------------------------------------------------------------------*

*--------------------------------------------------------------------------*
* (4) Standard error: placebo
*--------------------------------------------------------------------------*

*--------------------------------------------------------------------------*
* (5) Return output
*--------------------------------------------------------------------------*
ereturn clear
ereturn scalar tau = `tau_hat'
ereturn scalar N = `N'
ereturn scalar T = `T'
ereturn scalar N_treated = `n_treated'
ereturn scalar lambda_unit = `lu'
ereturn scalar lambda_time = `lt'

ereturn local cmd "trop"
ereturn local depvar "`1'"
ereturn local unitvar "`2'"
ereturn local timevar "`3'"
ereturn local treatvar "`4'"

di as txt _newline "Triply Robust Panel Estimator"
di as txt "--------------------------------"
di as txt "ATT estimate = " as result %9.6f e(tau)
end

*--------------------------------------------------------------------------*
* (6) Graphing (Currently no graphs required but can leave this space for later)
*--------------------------------------------------------------------------*

*------------------------------------------------------------------------------*
* (7) Stata Subroutines
*------------------------------------------------------------------------------*

*------------------------------------------------------------------------------*
* Mata functions
*------------------------------------------------------------------------------*

cap mata: mata drop trop_cell_weights()
cap mata: mata drop trop_fit_wls()
cap mata: mata drop trop_svt()
cap mata: mata drop trop_fit_nuclear()

mata:
void trop_cell_weights(
    real matrix Y,
    real matrix W,
    real colvector treated_units,
    real scalar lambda_unit,
    real scalar lambda_time,
    real matrix delta,
    real colvector delta_unit,
    real rowvector delta_time
)
{
    real scalar N, T, T0, center, treated_periods
    real rowvector dist_time, avg_treated, col_treated
    real colvector dist_unit, A, B
    real matrix mask, sqdiff

    N = rows(Y)
    T = cols(Y)

    if (rows(treated_units) == 0) {
        errprintf("No treated units supplied to trop_cell_weights().\n")
        exit(459)
    }

    col_treated     = colsum(W) :> 0
    treated_periods = sum(col_treated)
    T0              = T - treated_periods

    // Time distance: decay away from the middle of the treated block.
    center    = T - treated_periods / 2
    dist_time = abs(((1..T) :- 1) :- center)

    // Pre-period mask: 1 in pre-treatment columns, 0 in last treated_periods.
    mask = J(N, T, 1)
    mask[., (T0 + 1)..T] = J(N, treated_periods, 0)

    // Reference trajectory: average over treated units.
    avg_treated = mean(Y[treated_units, .])

    // Unit distance: RMS pre-period distance to the treated reference path.
    sqdiff    = (J(N, 1, 1) * avg_treated - Y) :^ 2
    A         = rowsum(sqdiff :* mask)
    B         = rowsum(mask)
    dist_unit = sqrt(A :/ B)

    delta_unit = exp(-lambda_unit :* dist_unit)
    delta_time = exp(-lambda_time :* dist_time)

    delta = delta_unit * delta_time
}
end

mata:
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

    b = qrsolve(XtWX, XtWy)

    mu    = b[1]
    alpha = 0 \ b[2..N]
    beta  = (0, b[(N + 1)..(N + T - 1)]')
    tau   = b[p]
}
end

mata:
real matrix trop_svt(real matrix Z, real scalar thr)
{
    real matrix U, Vt
    real colvector s, s_thr
    real scalar r

    fullsvd(Z, U, s, Vt)

    r     = rows(s)
    s_thr = (s :- thr) :* ((s :- thr) :> 0)

    return( (U[., 1..r] :* s_thr') * Vt[1..r, .] )
}

mata:
void trop_fit_nuclear(
    real matrix Y,
    real matrix W,
    real matrix delta,
    real scalar lambda_nn,
    real scalar tol,
    real scalar max_iter,
    real scalar tau,
    real scalar mu,
    real colvector alpha,
    real rowvector beta,
    real matrix L,
    real scalar iters
)
{
    real scalar N, T, NT, p, k, delta_max, Lip, step, thr
    real scalar t_mom, t_new, change, normL
    real colvector w_vec, y_vec, b
    real matrix D_unit, D_time, X, XtWX
    real matrix Z, L_new, fe, R, grad

    N  = rows(Y)
    T  = cols(Y)
    NT = N * T

    D_unit = I(N)[., 2..N] # J(T, 1, 1)
    D_time = J(N, 1, 1) # I(T)[., 2..T]

    p = 1 + (N - 1) + (T - 1) + 1
    X = J(NT, 1, 1), D_unit, D_time, vec(W')

    w_vec = vec(delta')
    XtWX  = quadcross(X, w_vec, X)

    delta_max = max(delta)
    if (delta_max <= 0) {
        errprintf("All delta weights are zero.\n")
        exit(498)
    }

    Lip  = 2 * delta_max
    step = 1 / Lip
    thr  = lambda_nn / Lip

    L     = J(N, T, 0)
    Z     = L
    t_mom = 1

    for (k = 1; k <= max_iter; k++) {
        y_vec = vec((Y - Z)')

        b  = qrsolve(XtWX, quadcross(X, w_vec, y_vec))
        fe = rowshape(X * b, N)

        tau = b[p]

        R     = Y - fe
        grad  = -2 :* delta :* (R - Z)
        L_new = trop_svt(Z - step :* grad, thr)

        t_new = (1 + sqrt(1 + 4 * t_mom^2)) / 2
        Z     = L_new + ((t_mom - 1) / t_new) :* (L_new - L)

        normL  = sqrt(sum(L_new :^ 2))
        change = sqrt(sum((L_new - L) :^ 2))

        L     = L_new
        t_mom = t_new

        if (change <= tol * max((1, normL))) break
    }

    mu    = b[1]
    alpha = 0 \ b[2..N]
    beta  = (0, b[(N + 1)..(N + T - 1)]')
    iters = min((k, max_iter))
}
end
