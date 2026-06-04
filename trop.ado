*! trop: Triply Robust Panel Estimators
*! Version 0.x.x May 5, 2025 
*! Author: Clarke Damian, Justin Waddy
*! dclarke@fen.uchile.cl, j.waddy@exeter.ac.uk

/*
Versions
0.1.0 May 26, 2026: Initial release. Weighted TWFE + nuclear-norm low-rank
       (FISTA with adaptive restart); placebo cross-validation
       (cycle/joint, resample/k-fold); stratified block-bootstrap SE.
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
    nn_grid(string)
    kfold(integer 0)
    level(integer 95)
    returnweights
    generate(string)
    verbose
    ntrials(integer 0)
    ntreated(integer 0)
    cv_seed(integer 0)
    ]
    ;
#delimit cr


*------------------------------------------------------------------------------*
* (0) Error checks
*------------------------------------------------------------------------------*
tempvar touse
mark `touse' `if' `in'

local unitname : word 2 of `varlist'

local stringvar = 0
capture confirm string variable `unitname'
if !_rc {
    local stringvar = 1
    tempvar ID
    quietly egen `ID' = group(`unitname')
    local varlist `1' `ID' `3' `4'
}
tokenize `varlist'

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

//DC: This restore below seems like it serves no purpose but adds time.  Can we confirm?
//    The reason is that we "re-preserve" on line 179, but in between we do nothign that
//    requires original data.  Seems more efficient to just stay preserved and pull out
//    block on line 179
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

if "`vce'" == "" local vce "bootstrap"
if !inlist("`vce'", "bootstrap", "noinference") {
    di as error "vce() must be one of: bootstrap, noinference."
    exit 198
}

*------------------------------------------------------------------------------*
* (1) Set-up 
*------------------------------------------------------------------------------*
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
if `lambda_unit_set' local lu = `lambda_unit'
else                 local lu = 0

if `lambda_time_set' local lt = `lambda_time'
else                 local lt = 0

* lambda_nn: set-and-inf -> missing (.) ; set-and-finite -> value ; unset -> 0 placeholder
if `lambda_nn_set' & `lambda_nn_inf'  local lnn = .
else if `lambda_nn_set'               local lnn = `lambda_nn'
else                                  local lnn = .

*--- Default CV grids (overridable via options) ---*
local default_unit_grid "0 0.1 0.2 0.3 0.5 0.8 1.2 1.6 2"
local default_time_grid "0 0.025 0.05 0.1 0.2 0.35 0.5 0.75 1 2 4"
local default_nn_grid   "0.005 0.01 0.025 0.05 0.1 0.25 0.5 1 ."

if "`unit_grid'" == "" local unit_grid "`default_unit_grid'"
if "`time_grid'" == "" local time_grid "`default_time_grid'"
if "`nn_grid'"   == "" local nn_grid   "`default_nn_grid'"

*K for placebo folds (default 5, overridable via kfold())*
if `kfold' <= 0 local kfold = 5

*--- CV method/search parsing: cv("method search"), defaults resample cycle ---*
local cv_method "resample"
local cv_search "cycle"
if "`cv'" != "" {
    tokenize "`cv'"
    if "`1'" != "" local cv_method "`1'"
    if "`2'" != "" local cv_search "`2'"
    * restore positional macros clobbered by tokenize
    tokenize `varlist'
}
if !inlist("`cv_method'", "resample", "kfold") {
    di as error "cv(): first word (method) must be 'resample' or 'kfold'."
    exit 198
}
if !inlist("`cv_search'", "cycle", "joint") {
    di as error "cv(): second word (search) must be 'cycle' or 'joint'."
    exit 198
}

if `ntrials'  <= 0 local ntrials  = 200
if `ntreated' <= 0 local ntreated = 1

*Run CV if any lambda is unspecified*
local need_cv = (1 - `lambda_unit_set') + (1 - `lambda_time_set') + (1 - `lambda_nn_set')

if `need_cv' > 0 {
    * For JOINT search, honor any user-fixed lambda by collapsing its grid to
    * the single fixed value (joint has no per-coordinate tune flags).
    if "`cv_search'" == "joint" {
        if `lambda_unit_set' local unit_grid "`lu'"
        if `lambda_time_set' local time_grid "`lt'"
        if `lambda_nn_set' {
            if "`lnn'" == "." local nn_grid "."
            else              local nn_grid "`lnn'"
        }
    }

    if "`verbose'" != "" {
        di as txt "Cross-validating lambdas (`cv_method' `cv_search', seed `cv_seed')..."
    }

    * Seed in Stata-land (matches bootstrap pattern); Mata runiform shares the stream.
    set seed `cv_seed'

    mata: trop_cv_setup(Y, W)

    if "`verbose'" != "" {
        di as txt "  selected lambda_unit = `lu'"
        di as txt "  selected lambda_time = `lt'"
        di as txt "  selected lambda_nn   = `lnn'"
    }
}

*Estimate on the full panel with resolved lambdas*
mata: delta = .; delta_unit = .; delta_time = .
mata: trop_cell_weights(Y, W, strtoreal(tokens(st_local("treated_rows")))', `lu', `lt', delta, delta_unit, delta_time)

mata: tau_hat = .; mu_hat = .; alpha_hat = .; beta_hat = .; L_hat = .; iters = .
if "`lnn'" == "." {
    mata: trop_fit_wls(Y, W, delta, tau_hat, mu_hat, alpha_hat, beta_hat)
}
else {
    mata: trop_fit_nuclear(Y, W, delta, `lnn', 1e-10, 5000, tau_hat, mu_hat, alpha_hat, beta_hat, L_hat, iters)
}
mata: st_local("tau_hat", strofreal(tau_hat))


*--------------------------------------------------------------------------*
* (3) Standard error: bootstrap
*--------------------------------------------------------------------------*
local se = .
if "`vce'" != "noinference" {
    if `seed' != 0  set seed `seed'
    if "`verbose'" != "" di as txt "Bootstrapping SE (`reps' reps)..."

    if "`lnn'" == "." local lnn_arg = .
    else              local lnn_arg = `lnn'

    mata: st_numscalar("r(se)", trop_bootstrap_se(Y, W, `lu', `lt', `lnn_arg', 1e-10, 5000, `reps'))
    local se = r(se)
}

*--------------------------------------------------------------------------*
* (4) Return output
*--------------------------------------------------------------------------*
ereturn clear
ereturn scalar tau = `tau_hat'
ereturn scalar N = `N'
ereturn scalar T = `T'
ereturn scalar N_treated = `n_treated'
ereturn scalar lambda_unit = `lu'
ereturn scalar lambda_time = `lt'
if "`lnn'" == "." ereturn scalar lambda_nn = .
else              ereturn scalar lambda_nn = `lnn'

ereturn local cmd "trop"
ereturn local depvar "`1'"
ereturn local unitvar "`unitname'"
ereturn local timevar "`3'"
ereturn local treatvar "`4'"
ereturn local vce "`vce'"

* Inference scalars (only if we bootstrapped)
local have_se = ("`se'" != "." & "`se'" != "")
if `have_se' {
    ereturn scalar se = `se'
    local crit = invnormal(1 - (1 - `level'/100)/2)
    ereturn scalar ci_lower = `tau_hat' - `crit' * `se'
    ereturn scalar ci_upper = `tau_hat' + `crit' * `se'
    ereturn scalar z = `tau_hat' / `se'
    ereturn scalar pvalue = 2 * normal(-abs(`tau_hat' / `se'))
    ereturn scalar level = `level'
}

* lambda_nn display string
if "`lnn'" == "." local lnn_disp "inf"
else              local lnn_disp "`lnn'"

*--- Output table ---*
di ""
di as text "{hline 13}{c TT}{hline 50}"
di as text %12s "TROP" " {c |}  Triply Robust Panel estimator"
di as text "{hline 13}{c +}{hline 50}"
di as text %12s "ATT" " {c |}  " as result %10.5f e(tau)
if `have_se' {
    local cilab = "`level'% CI"
    di as text %12s "Std. err." " {c |}  " as result %10.5f e(se)
    di as text %12s "z" " {c |}  " as result %10.4f e(z)
    di as text %12s "P>|z|" " {c |}  " as result %10.4f e(pvalue)
    di as text %12s "`cilab'" " {c |}  " as result %9.5f e(ci_lower) "  " %9.5f e(ci_upper)
}
else {
    di as text %12s " " " {c |}  " as text "(no inference; vce(noinference))"
}
di as text "{hline 13}{c +}{hline 50}"
di as text %12s "N units" " {c |}  " as result %10.0f e(N)
di as text %12s "T periods" " {c |}  " as result %10.0f e(T)
di as text %12s "N treated" " {c |}  " as result %10.0f e(N_treated)
di as text "{hline 13}{c +}{hline 50}"
di as text %12s "lambda_unit" " {c |}  " as result %10.4f `lu'
di as text %12s "lambda_time" " {c |}  " as result %10.4f `lt'
di as text %12s "lambda_nn" " {c |}  " as result %10s "`lnn_disp'"
if `need_cv' > 0 {
    di as text %12s " " " {c |}  " as text "(selected by `kfold'-fold CV)"
}
di as text "{hline 13}{c BT}{hline 50}"
end

*------------------------------------------------------------------------------*
* Mata functions
*------------------------------------------------------------------------------*

cap mata: mata drop trop_cell_weights()
cap mata: mata drop trop_fit_wls()
cap mata: mata drop trop_svt()
cap mata: mata drop trop_fit_nuclear()
cap mata: mata drop trop_placebo_rmse()
cap mata: mata drop trop_cv_single()
cap mata: mata drop trop_cv_cycle()
cap mata: mata drop trop_cv_setup()
cap mata: mata drop trop_kfold_sets()
cap mata: mata drop trop_resample_sets()
cap mata: mata drop trop_cv_joint()

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
real matrix trop_svt(real matrix Z, real scalar thr, real scalar nucnorm)
{
    real matrix U, Vt
    real colvector s, s_thr
    real scalar r

    fullsvd(Z, U, s, Vt)

    r       = rows(s)
    s_thr   = (s :- thr) :* ((s :- thr) :> 0)
    nucnorm = sum(s_thr)                 // nuclear norm of the output = sum of thresholded singular values
    return( (U[., 1..r] :* s_thr') * Vt[1..r, .] )
}
end

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
    real scalar t_mom, t_new
    real scalar tol_tau, tol_L, tol_obj
    real scalar tau_new, tau_old, obj_new, obj_old, nucnorm, datafit
    real scalar tau_change, L_change, obj_change, normL, restart_grad
    real colvector w_vec, y_vec, b
    real matrix D_unit, D_time, X, XtWX
    real matrix Z, L_new, fe, R, grad

    N  = rows(Y); T = cols(Y); NT = N * T

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

    // Fixed step from a global Lipschitz bound on the smooth weighted loss.
    Lip  = 2 * delta_max
    step = 1 / Lip
    thr  = lambda_nn / Lip

    // Tolerances derived from the single `tol`. The treatment estimate (tau) is
    // the primary convergence criterion; the low-rank component and penalized
    // objective are looser "still-moving?" safeguards, so we do not stop while
    // L is still rearranging the factor structure (which can leave tau briefly
    // stalled and then move again).
    tol_tau = tol
    tol_L   = tol;  if (tol_L   < 1e-8) tol_L   = 1e-8
    tol_obj = tol;  if (tol_obj < 1e-8) tol_obj = 1e-8

    L       = J(N, T, 0)
    Z       = L
    t_mom   = 1
    tau_old = .         // missing -> first-iteration changes are non-binding
    obj_old = .
    nucnorm = 0

    for (k = 1; k <= max_iter; k++) {

        // (a) Fixed-effect step: exact weighted LS on (Y - Z), at extrapolation Z.
        y_vec   = vec((Y - Z)')
        b       = qrsolve(XtWX, quadcross(X, w_vec, y_vec))
        fe      = rowshape(X * b, N)
        tau     = b[p]
        tau_new = tau

        // (b) Low-rank proximal-gradient (FISTA) step; gradient taken at Z.
        R     = Y - fe
        grad  = -2 :* delta :* (R - Z)
        L_new = trop_svt(Z - step :* grad, thr, nucnorm)

        // (c) Penalized objective at the current iterate (fe, L_new). Data fit
        //     reuses R; nuclear norm comes back from trop_svt
        datafit = sum(delta :* ((R - L_new):^2))
        obj_new = datafit + lambda_nn * nucnorm

        // (d) Adaptive restart (O'Donoghue & Candes 2015). Reset momentum if the
        //     gradient-mapping direction opposes progress (gradient scheme) OR
        //     the penalized objective increased (function scheme). Either flags
        //     an overshoot; resetting preserves the O(1/k^2) guarantee while
        //     suppressing the oscillation that hurts ill-conditioned weights.
        restart_grad = sum((Z - L_new) :* (L_new - L))
        if (restart_grad > 0 | obj_new > obj_old + 1e-10 * max((1, abs(obj_old)))) {
            t_new = 1
            Z     = L_new
        }
        else {
            t_new = (1 + sqrt(1 + 4 * t_mom^2)) / 2
            Z     = L_new + ((t_mom - 1) / t_new) :* (L_new - L)
        }

        // (e) Relative changes for the stopping rule.
        normL      = sqrt(sum(L_new:^2))
        tau_change = abs(tau_new - tau_old) / max((1, abs(tau_new)))
        L_change   = sqrt(sum((L_new - L):^2)) / max((1, normL))
        obj_change = abs(obj_new - obj_old) / max((1, abs(obj_new)))

        // Advance state.
        L       = L_new
        t_mom   = t_new
        tau_old = tau_new
        obj_old = obj_new

        // (f) Converged when the treatment estimate is stable AND the low-rank
        //     component and objective are no longer materially moving.
        if (tau_change < tol_tau & L_change < tol_L & obj_change < tol_obj) break
    }

    mu    = b[1]
    alpha = 0 \ b[2..N]
    beta  = (0, b[(N + 1)..(N + T - 1)]')
    iters = min((k, max_iter))
}
end

// One TROP placebo estimate given a placebo W on the control panel.
// Dispatches WLS (lambda_nn = .) vs nuclear, exactly as section (2) does.
mata:
real scalar trop_placebo_tau(
    real matrix Yc, real matrix Wp,
    real scalar lu, real scalar lt, real scalar lnn,
    real scalar tol, real scalar max_iter)
{
    real matrix delta, L_hat
    real colvector delta_unit, alpha_hat, treated_units
    real rowvector delta_time, beta_hat
    real scalar tau_hat, mu_hat, iters

    delta = .; delta_unit = .; delta_time = .
    treated_units = selectindex(rowsum(Wp) :> 0)
    trop_cell_weights(Yc, Wp, treated_units, lu, lt, delta, delta_unit, delta_time)

    tau_hat = .; mu_hat = .; alpha_hat = .; beta_hat = .; L_hat = .; iters = .
    if (lnn >= .) {     // missing encodes lambda_nn = infinity (WLS path)
        trop_fit_wls(Yc, Wp, delta, tau_hat, mu_hat, alpha_hat, beta_hat)
    }
    else {
        trop_fit_nuclear(Yc, Wp, delta, lnn, tol, max_iter,
                         tau_hat, mu_hat, alpha_hat, beta_hat, L_hat, iters)
    }
    return(tau_hat)
}
end


// Placebo RMSE for one lambda triplet, given pre-built placebo SETS.
// `sets` is a pointer rowvector; each *sets[j] is a colvector of control-panel
// row indices to treat (in the last treated_periods columns) for placebo j.
// Identical for k-fold sets (each unit in exactly one set) and resample sets
// (n_trials independent draws). Returns sqrt(mean(tau^2)) over finite estimates.
mata:
real scalar trop_placebo_rmse(
    real matrix Yc, real scalar treated_periods,
    pointer(real colvector) rowvector sets,
    real scalar lu, real scalar lt, real scalar lnn,
    real scalar tol, real scalar max_iter)
{
    real scalar Nc, T, J_sets, j, tau, nfin, ssq
    real colvector set_rows
    real matrix Wp

    Nc = rows(Yc)
    T  = cols(Yc)
    J_sets = cols(sets)

    nfin = 0
    ssq  = 0
    for (j = 1; j <= J_sets; j++) {
        set_rows = *sets[j]                         
        Wp = J(Nc, T, 0)
        Wp[set_rows, (T - treated_periods + 1)..T] = J(rows(set_rows), treated_periods, 1)
        tau = trop_placebo_tau(Yc, Wp, lu, lt, lnn, tol, max_iter)
        if (!missing(tau)) {
            ssq  = ssq + tau^2
            nfin = nfin + 1
        }
    }
    if (nfin == 0) return(.)
    return(sqrt(ssq / nfin))
}
end

mata:
// cv_single: tune one lambda over a grid (placebo CV), others fixed.
//   which_lambda : 1 = unit, 2 = time, 3 = nn
//   fixed1, fixed2: the other two, in canonical order
// Sets are PASSED IN (built once by the caller) so every lambda on the grid —
// and every coordinate/cycle in cv_cycle — is scored on the SAME placebo folds.
// This matches the Python reference (draw once, reuse) and makes cross-lambda
// RMSE comparisons apples-to-apples.
real scalar trop_cv_single(
    real matrix Yc,
    real scalar treated_periods,
    pointer(real colvector) rowvector sets,
    real rowvector grid,
    real scalar which_lambda,
    real scalar fixed1,
    real scalar fixed2,
    real scalar tol,
    real scalar max_iter,
    real rowvector scores
)
{
    real scalar ng, g, lu, lt, lnn
    real scalar best_idx, best_score

    ng = cols(grid)
    scores = J(1, ng, .)

    for (g = 1; g <= ng; g++) {
        if (which_lambda == 1) {
            lu = grid[g]; lt = fixed1; lnn = fixed2
        }
        else if (which_lambda == 2) {
            lu = fixed1; lt = grid[g]; lnn = fixed2
        }
        else if (which_lambda == 3) {
            lu = fixed1; lt = fixed2; lnn = grid[g]
        }
        else {
            errprintf("trop_cv_single(): which_lambda must be 1, 2, or 3.\n")
            exit(198)
        }

        scores[g] = trop_placebo_rmse(
            Yc, treated_periods, sets, lu, lt, lnn, tol, max_iter)
    }

    best_idx = .; best_score = .
    for (g = 1; g <= ng; g++) {
        if (scores[g] < .) {
            if (best_score >= . | scores[g] < best_score) {
                best_score = scores[g]; best_idx = g
            }
        }
    }
    if (best_idx >= .) {
        errprintf("trop_cv_single(): all CV scores are missing.\n")
        exit(498)
    }
    return(grid[best_idx])
}
end

mata:
// cv_cycle with per-coordinate tune flags. Sets are PASSED IN (built once by
// cv_setup, by the chosen sampling method) and reused across all coordinates
// and cycles. A coordinate updates each cycle only if its flag is nonzero;
// fixed coordinates stay at their passed-in value.
real scalar trop_cv_cycle(
    real matrix Yc, real scalar treated_periods,
    pointer(real colvector) rowvector sets,
    real rowvector unit_grid, real rowvector time_grid, real rowvector nn_grid,
    real scalar tune_unit, real scalar tune_time, real scalar tune_nn,
    real scalar tol, real scalar max_iter, real scalar max_cycles,
    real scalar lambda_unit, real scalar lambda_time, real scalar lambda_nn,
    real scalar cycles_used)
{
    real scalar c, lu_old, lt_old, lnn_old
    real rowvector scores, fin

    // Initialize tuned coordinates to grid finite-mean; leave fixed ones as passed in.
    if (tune_unit) {
        fin = select(unit_grid, unit_grid :< .)
        lambda_unit = mean(fin')
    }
    if (tune_time) {
        fin = select(time_grid, time_grid :< .)
        lambda_time = mean(fin')
    }
    if (tune_nn) {
        fin = select(nn_grid, nn_grid :< .)
        lambda_nn = mean(fin')
    }

    scores = .
    for (c = 1; c <= max_cycles; c++) {
        lu_old = lambda_unit; lt_old = lambda_time; lnn_old = lambda_nn

        if (tune_unit)
            lambda_unit = trop_cv_single(Yc, treated_periods, sets, unit_grid, 1,
                                         lambda_time, lambda_nn, tol, max_iter, scores)
        if (tune_time)
            lambda_time = trop_cv_single(Yc, treated_periods, sets, time_grid, 2,
                                         lambda_unit, lambda_nn, tol, max_iter, scores)
        if (tune_nn)
            lambda_nn = trop_cv_single(Yc, treated_periods, sets, nn_grid, 3,
                                       lambda_unit, lambda_time, tol, max_iter, scores)

        if (lambda_unit == lu_old & lambda_time == lt_old & lambda_nn == lnn_old) {
            cycles_used = c
            return(lambda_unit)
        }
    }
    cycles_used = max_cycles
    return(lambda_unit)
}
end

mata:
void trop_cv_setup(real matrix Yfull, real matrix Wfull)
{
    real matrix Yc
    real rowvector unit_grid, time_grid, nn_grid
    real colvector control_rows
    real scalar treated_periods, K, tol, max_iter, max_cycles, Nc
    real scalar tune_unit, tune_time, tune_nn
    real scalar lambda_unit, lambda_time, lambda_nn, cycles_used
    real scalar ntrials, ntreated, best_rmse, n_eval
    string scalar method, search
    pointer(real colvector) rowvector sets

    // Control-only panel: rows never treated.
    control_rows = selectindex(rowsum(Wfull) :== 0)
    if (rows(control_rows) < 2) {
        errprintf("Need at least 2 never-treated units for cross-validation.\n")
        exit(459)
    }
    Yc = Yfull[control_rows, .]
    Nc = rows(Yc)

    treated_periods = sum(colsum(Wfull) :> 0)

    unit_grid = strtoreal(tokens(st_local("unit_grid")))
    time_grid = strtoreal(tokens(st_local("time_grid")))
    nn_grid   = strtoreal(tokens(st_local("nn_grid")))

    K          = strtoreal(st_local("kfold"))
    tol        = 1e-10
    max_iter   = 5000
    max_cycles = 50

    tune_unit = 1 - strtoreal(st_local("lambda_unit_set"))
    tune_time = 1 - strtoreal(st_local("lambda_time_set"))
    tune_nn   = 1 - strtoreal(st_local("lambda_nn_set"))

    lambda_unit = strtoreal(st_local("lu"))
    lambda_time = strtoreal(st_local("lt"))
    lambda_nn   = strtoreal(st_local("lnn"))

    // CV method/search and resample params (defaults set in the command body).
    method   = st_local("cv_method")
    search   = st_local("cv_search")
    ntrials  = strtoreal(st_local("ntrials"))
    ntreated = strtoreal(st_local("ntreated"))

    // Build placebo sets ONCE, by the chosen sampling method. The command body
    // has already `set seed`'d, so these draws are reproducible.
    if (method == "resample") {
        sets = trop_resample_sets(Nc, ntrials, ntreated)
    }
    else {   // "kfold"
        sets = trop_kfold_sets(Nc, K)
    }

    cycles_used = .
    if (search == "joint") {
        // Exhaustive grid. Joint always tunes all three dimensions over the
        // supplied grids; fixed coordinates are honored by passing a 1-element
        // grid (handled in the command body when a lambda is user-set).
        best_rmse = .; n_eval = .
        trop_cv_joint(Yc, treated_periods, sets,
                      unit_grid, time_grid, nn_grid,
                      tol, max_iter,
                      lambda_unit, lambda_time, lambda_nn, best_rmse, n_eval)
        st_local("cv_cycles", "0")   // joint is not iterative
    }
    else {   // "cycle"
        (void) trop_cv_cycle(Yc, treated_periods, sets,
                             unit_grid, time_grid, nn_grid,
                             tune_unit, tune_time, tune_nn,
                             tol, max_iter, max_cycles,
                             lambda_unit, lambda_time, lambda_nn, cycles_used)
        st_local("cv_cycles", strofreal(cycles_used))
    }

    st_local("lu", strofreal(lambda_unit))
    st_local("lt", strofreal(lambda_time))
    if (lambda_nn >= .) st_local("lnn", ".")
    else                st_local("lnn", strofreal(lambda_nn))
}
end

cap mata: mata drop trop_bootstrap_se()

mata:

// Stratified bootstrap SE (Algorithm 3). Resample N1 treated rows and N0
// control rows with replacement, rebuild the panel, recompute tau at FIXED
// lambdas, repeat reps times. SE = sqrt((B-1)/B) * sd(tau_b).
//
// Seeding: caller does `set seed` in Stata before this is invoked; Mata's
// runiform() shares Stata's RNG stream, so the draws are reproducible.
real scalar trop_bootstrap_se(
    real matrix Y, real matrix W,
    real scalar lu, real scalar lt, real scalar lnn,
    real scalar tol, real scalar max_iter, real scalar reps)
{
    real scalar b, tau_b, mu_b, iters
    real scalar B_eff, tbar, ssq
    real colvector ctrl_rows, trt_rows, n0, n1
    real colvector samp_ctrl, samp_trt, samp_rows
    real colvector alpha_b, delta_unit, treated_units
    real rowvector beta_b, delta_time
    real matrix Yb, Wb, delta, L_b
    real colvector tau_store

    // Row partition by ever-treated status.
    trt_rows  = selectindex(rowsum(W) :> 0)
    ctrl_rows = selectindex(rowsum(W) :== 0)
    n1 = rows(trt_rows)
    n0 = rows(ctrl_rows)

    if (n1 == 0 | n0 == 0) {
        errprintf("Bootstrap requires both treated and control units.\n")
        return(.)
    }

    tau_store = J(reps, 1, .)

    for (b = 1; b <= reps; b++) {
        // Resample WITH replacement within each stratum.
        samp_trt  = trt_rows[ ceil(n1 :* runiform(n1, 1)) ]
        samp_ctrl = ctrl_rows[ ceil(n0 :* runiform(n0, 1)) ]
        samp_rows = samp_trt \ samp_ctrl

        Yb = Y[samp_rows, .]
        Wb = W[samp_rows, .]

        // Recompute weights and estimate at FIXED lambdas.
        delta = .; delta_unit = .; delta_time = .
        treated_units = selectindex(rowsum(Wb) :> 0)
        if (rows(treated_units) == 0) continue   // degenerate draw, skip

        trop_cell_weights(Yb, Wb, treated_units, lu, lt, delta, delta_unit, delta_time)

        tau_b = .; mu_b = .; alpha_b = .; beta_b = .; L_b = .; iters = .
        if (lnn >= .) {
            trop_fit_wls(Yb, Wb, delta, tau_b, mu_b, alpha_b, beta_b)
        }
        else {
            trop_fit_nuclear(Yb, Wb, delta, lnn, tol, max_iter,
                             tau_b, mu_b, alpha_b, beta_b, L_b, iters)
        }
        tau_store[b] = tau_b
    }

    // SE over finite draws: sqrt((B-1)/B) * sample sd.
    tau_store = select(tau_store, tau_store :< .)
    B_eff = rows(tau_store)
    if (B_eff < 2) {
        errprintf("Too few finite bootstrap replications for an SE.\n")
        return(.)
    }
    tbar = mean(tau_store)
    ssq  = sum((tau_store :- tbar):^2)
    // sample variance uses (B_eff - 1); the (B-1)/B factor matches the
    // convention in the transcript / sdid.
    return( sqrt((B_eff - 1) / B_eff) * sqrt(ssq / (B_eff - 1)) )
}

end

mata:

// Seeding: caller does `set seed` in Stata first; Mata runiform() shares the
// stream, so the permutation is deterministic given that seed.
pointer(real colvector) rowvector trop_kfold_sets(real scalar Nc, real scalar K)
{
    real colvector perm
    real scalar base, rem, f, start, len
    pointer(real colvector) rowvector sets

    perm = order(runiform(Nc, 1), 1)        // seeded uniform random permutation
    sets = J(1, K, NULL)
    base = floor(Nc / K)
    rem  = Nc - base * K
    start = 1
    for (f = 1; f <= K; f++) {
        len = base + (f <= rem)             // first `rem` folds get the extra unit
        sets[f] = &(perm[start..(start + len - 1)])
        start = start + len
    }
    return(sets)
}
end

mata:
pointer(real colvector) rowvector trop_resample_sets(
    real scalar Nc, real scalar n_trials, real scalar n_treated)
{
    real scalar i
    real colvector perm
    pointer(real colvector) rowvector sets

    if (n_treated > Nc) {
        errprintf("trop_resample_sets(): n_treated (%g) exceeds Nc (%g).\n",
                  n_treated, Nc)
        exit(198)
    }

    sets = J(1, n_trials, NULL)
    for (i = 1; i <= n_trials; i++) {
        perm    = order(runiform(Nc, 1), 1)   // fresh seeded permutation each trial
        sets[i] = &(perm[1..n_treated])       // first n_treated = w/o-replacement draw
    }
    return(sets)
}
end

mata:
void trop_cv_joint(
    real matrix Yc, real scalar treated_periods,
    pointer(real colvector) rowvector sets,
    real rowvector unit_grid, real rowvector time_grid, real rowvector nn_grid,
    real scalar tol, real scalar max_iter,
    real scalar lambda_unit, real scalar lambda_time, real scalar lambda_nn,
    real scalar best_rmse, real scalar n_eval)
{
    real scalar nu, nt, nn, ii, jj, kk, score
    real scalar best_u, best_t, best_n

    nu = cols(unit_grid)
    nt = cols(time_grid)
    nn = cols(nn_grid)

    best_rmse = .        
    best_u = .; best_t = .; best_n = .
    n_eval = 0

    for (ii = 1; ii <= nu; ii++) {
        for (jj = 1; jj <= nt; jj++) {
            for (kk = 1; kk <= nn; kk++) {

                score = trop_placebo_rmse(
                    Yc, treated_periods, sets,
                    unit_grid[ii], time_grid[jj], nn_grid[kk],
                    tol, max_iter)
                n_eval = n_eval + 1

                if (score < .) {
                    if (best_rmse >= . | score < best_rmse) {
                        best_rmse = score
                        best_u    = unit_grid[ii]
                        best_t    = time_grid[jj]
                        best_n    = nn_grid[kk]
                    }
                }
            }
        }
    }

    if (best_rmse >= .) {
        errprintf("trop_cv_joint(): all CV scores are missing.\n")
        exit(498)
    }

    lambda_unit = best_u
    lambda_time = best_t
    lambda_nn   = best_n     // may be missing (.) -> inf/WLS, correct
}
end
