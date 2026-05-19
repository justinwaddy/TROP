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
    treated_periods(integer 10)
    solver(string)
    vce(string)
    reps(integer 200)
    seed(integer 0)
    cv(string)
    unit_grid(numlist)
    time_grid(numlist)
    nn_grid(numlist)
    n_treated_units(integer 1)
    kfold(integer 0)
    n_jobs(integer 1)
    level(integer 95)
    graph
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
local N = r(N_g)
local T = r(Tmax)

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
*Core function will be in this section

*First subsection will be cross-validation code to tune parameters for lambda_unit, lambda_time and lambda_nn if they are not provided.
*If these are provided, then next subsection determines the ATT

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

*--------------------------------------------------------------------------*
* (6) Graphing (Currently no graphs required but can leave this space for later)
*--------------------------------------------------------------------------*

*------------------------------------------------------------------------------*
* (7) Stata Subroutines
*------------------------------------------------------------------------------*

*------------------------------------------------------------------------------*
* (8) Mata functions
*------------------------------------------------------------------------------*
cap mata: mata drop trop_cell_weights()

mata:
void trop_cell_weights(Y, treated_units, lambda_unit, lambda_time,
                         treated_periods, delta, delta_unit, delta_time)
{
    real scalar N, T, center
    real rowvector dist_time, avg_treated
    real colvector dist_unit, A, B
    real matrix mask, sqdiff

    N = rows(Y)
    T = cols(Y)

    // Time distances (matches Python's 0-based np.arange convention)
    center    = T - treated_periods / 2
    dist_time = abs(((1..T) :- 1) :- center)

    // Pre-period mask: 1 in pre, 0 in last treated_periods columns
    mask = J(N, T, 1)
    mask[., (T - treated_periods + 1)..T] = J(N, treated_periods, 0)

    // Reference trajectory: column mean over treated-unit rows
    avg_treated = mean(Y[treated_units, .])

    // RMS pre-period distance from each unit to the reference
    sqdiff    = (J(N, 1, 1) * avg_treated - Y) :^ 2
    A         = rowsum(sqdiff :* mask)
    B         = rowsum(mask)
    dist_unit = sqrt(A :/ B)

    delta_unit = exp(-lambda_unit :* dist_unit)
    delta_time = exp(-lambda_time :* dist_time)
    delta      = delta_unit * delta_time
}
end
