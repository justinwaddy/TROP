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
version 18.0

#delimit ;
syntax varlist(min=4 max=4) [if] [in],
    [
    lambda_unit(real 0)
    lambda_time(real 0)
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
*Write this after building the core functions

*------------------------------------------------------------------------------*
* (1) Set-up 
*------------------------------------------------------------------------------*
*Prepare dataset before treatment. See example prep: https://github.com/ostasovskyi/TROP-Estimator/blob/main/notebooks/tutorial.ipynb?short_path=b4fd5c7

if (length("`if'")+length("`in'")>0) {
    preserve
    qui keep if `touse'
}

*------------------------------------------------------------------------------*
* (2) Calculate ATT (jackknife will be default)
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