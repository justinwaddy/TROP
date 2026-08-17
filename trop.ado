*! trop: Triply Robust Panel Estimators
*! Version 0.2.8 August 17, 2026
*! Author: Clarke Damian, Justin Waddy
*! dclarke@fen.uchile.cl, j.waddy@exeter.ac.uk

/*
Versions
0.2.8 August 13, 2026: Adaptive CV now begins with a sweep of the full lambda grids (the defaults are the same
      grids used by cycle and joint). Verbose heatmap now has a Legend, which maps each symbol to its RMSE range.
0.2.7 August 7, 2026: Added adaptive CV, vce(jackknife) and progress bars. Verbose shows ASCII visual of adaptive CV.
0.2.6 August 4, 2026: Pooled (group(time)) fixed to correctly mask other treated units/current time period for unit distances.
      Added pooled_treat_distance(time() unit()) option for group(time): allows you to set how distances to own treated group
      are calculated. Default is empty, or ".", which uses the usual distances (i.e. from  midpoint for time weights or from
      average of outcomes at each time period across the i units in a group for unit weights). Setting a number, such as zero, 
      overrides the distances for treated units.
0.2.5 July 6, 2026: Added 'detail' option which provides further details on heterogeneous treatment effects and the units in a given cohort.
0.2.4 July 4, 2026: Resample and K-fold now use the pattern of treatment from treated units instead of entire block. Pooled time weights 
      now center on the block midpoint (min+max)/2, aligning with dist(s,t)=|t-s| (previously half-period late). Single treated 
      periods with the same lambdas now result in the same estimates under group(time) and group(cell).
0.2.3 June 29, 2026: Default settings changed from group(cell) to group(time), and no inference is computed by default. 
0.2.2 June 11, 2026: Fixed LOOCV to be per-cell (fixes group(cell)). Added cells() suboption which specifies
      how many cells you want to randomly sample using LOOCV for large panels. Added additional validation checks on CV options
      depending on grouping by cell or time.
0.2.1 June 4, 2026: Standardized the outcome ((Y - overall mean)/overall SD)
      before fitting; tau/SE/CI are mapped back to the raw scale. Rescale makes lambda_unit/lambda_time/
      lambda_nn grids standardized i.e. returned lambdas are the standardized scale.
0.2.0 June 3, 2026: Optimised FISTA, added warm-start nuclear path, and optimised CV algorithms.
      Added generalized treatment patterns and bootstrapping. FISTA optimises when tau converges rather than
      L. 
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
    group(string)
    lambda_unit(string)
    lambda_time(string)
    lambda_nn(string)
    cv(string)
    vce(string)
    pooled_treat_distance(string)
    level(integer 95)
    verbose
    detail
    ]
    ;
#delimit cr


*------------------------------------------------------------------------------*
* (0) Error checks
*------------------------------------------------------------------------------*

*--- group() ---*
if "`group'" == "" local group "time"
if !inlist("`group'", "cell", "time") {
    di as error "group() must be 'cell' or 'time'."
    exit 198
}

if "`group'" == "time" local cv_method "resample"
else                   local cv_method "loocv"
local cv_search "cycle"
local cv_seed   0
local ntrials   200
local kfold     5
local ncells    0
local npoints     10
local nexpansions 6
local unit_grid ""
local time_grid ""
local nn_grid   ""
if `"`cv'"' != "" {
    local _cc = strpos(`"`cv'"', ",")
    if `_cc' == 0 local _cvmain = trim(`"`cv'"')
    else          local _cvmain = trim(substr(`"`cv'"', 1, `_cc'-1))
    local _m : word 1 of `_cvmain'
    local _s : word 2 of `_cvmain'
    if "`_m'" != "" local cv_method "`_m'"
    if "`_s'" != "" local cv_search "`_s'"
    if `_cc' != 0 {
        local _cvsub = substr(`"`cv'"', `_cc'+1, .)
        foreach _k in trials folds seed cells points expansions unit_grid time_grid nn_grid {
            local _p = strpos(`"`_cvsub'"', "`_k'(")
            if `_p' > 0 {
                local _rest = substr(`"`_cvsub'"', `_p' + length("`_k'("), .)
                local _q = strpos(`"`_rest'"', ")")
                local _`_k' = substr(`"`_rest'"', 1, `_q'-1)
            }
        }
        if "`_trials'"      != "" local ntrials     = `_trials'
        if "`_folds'"       != "" local kfold       = `_folds'
        if "`_seed'"        != "" local cv_seed     = `_seed'
        if "`_cells'"       != "" local ncells      = `_cells'
        if "`_points'"      != "" local npoints     = `_points'
        if "`_expansions'"  != "" local nexpansions = `_expansions'
        if `"`_unit_grid'"' != "" local unit_grid `"`_unit_grid'"'
        if `"`_time_grid'"' != "" local time_grid `"`_time_grid'"'
        if `"`_nn_grid'"'   != "" local nn_grid   `"`_nn_grid'"'
    }
}
if !inlist("`cv_method'", "loocv", "resample", "kfold") {
    di as error "cv(): method must be 'loocv', 'resample', or 'kfold'."
    exit 198
}
if !inlist("`cv_search'", "cycle", "joint", "adaptive") {
    di as error "cv(): search must be 'cycle', 'joint', or 'adaptive'."
    exit 198
}
if "`cv_search'" == "adaptive" {
    if `npoints' < 2 {
        di as error "cv(): points() must be at least 2."
        exit 198
    }
    if `nexpansions' < 0 {
        di as error "cv(): expansions() must be nonnegative."
        exit 198
    }
}

local reps 200
local seed 0
if `"`vce'"' != "" {
    local _vc = strpos(`"`vce'"', ",")
    if `_vc' == 0 local vce = trim(`"`vce'"')
    else {
        local _vct = trim(substr(`"`vce'"', 1, `_vc'-1))
        local _vcs = substr(`"`vce'"', `_vc'+1, .)
        foreach _k in reps seed {
            local _p = strpos(`"`_vcs'"', "`_k'(")
            if `_p' > 0 {
                local _rest = substr(`"`_vcs'"', `_p' + length("`_k'("), .)
                local _q = strpos(`"`_rest'"', ")")
                local _v`_k' = substr(`"`_rest'"', 1, `_q'-1)
            }
        }
        if "`_vreps'" != "" local reps = `_vreps'
        if "`_vseed'" != "" local seed = `_vseed'
        local vce "`_vct'"
    }
}
if "`vce'" == "" local vce "noinference"
if !inlist("`vce'", "bootstrap", "jackknife", "noinference") {
    di as error "vce() must be 'bootstrap', 'jackknife', or 'noinference'."
    exit 198
}
if "`vce'" == "jackknife" {
    if "`_vreps'" != "" di as txt "note: reps() ignored under vce(jackknife) (one replication per unit)."
    if "`_vseed'" != "" di as txt "note: seed() ignored under vce(jackknife) (deterministic)."
}

local ptd_unit_val "."
local ptd_time_val "."
if `"`pooled_treat_distance'"' != "" {
    local _ptd_time ""
    local _ptd_unit ""
    foreach _k in time unit {
        local _p = strpos(`"`pooled_treat_distance'"', "`_k'(")
        if `_p' > 0 {
            local _rest = substr(`"`pooled_treat_distance'"', `_p' + length("`_k'("), .)
            local _q = strpos(`"`_rest'"', ")")
            local _ptd_`_k' = trim(substr(`"`_rest'"', 1, `_q'-1))
        }
    }
    foreach _k in time unit {
        if "`_ptd_`_k''" != "" & "`_ptd_`_k''" != "." {
            cap confirm number `_ptd_`_k''
            if _rc {
                di as error "pooled_treat_distance(): `_k'() must be a number or '.'."
                exit 198
            }
        }
    }
    if "`group'" != "time" {
        di as txt "note: pooled_treat_distance() ignored under group(cell)."
    }
    else {
        if "`_ptd_unit'" != "" & "`_ptd_unit'" != "." local ptd_unit_val = `_ptd_unit'
        if "`_ptd_time'" != "" & "`_ptd_time'" != "." local ptd_time_val = `_ptd_time'
    }
}

tempvar touse
mark `touse' `if' `in'

local unitname : word 2 of `varlist'

capture confirm string variable `unitname'
if !_rc {
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
qui count if `4'==1
local ntrcells = r(N)
qui sum `4'
if (r(min)==0 & r(max)==0) {
    di as error "All units are controls."
    exit 459
}
qui sum `4'
if (r(min)==1 & r(max)==1) {
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
if (length("`if'")+length("`in'")>0) {
    preserve
    qui keep if `touse'
}

qui sort `2' `3'

*Wide panels
mata: Y = rowshape(st_data(., "`1'"), `N')
mata: W = rowshape(st_data(., "`4'"), `N')

*Pooled treated-block distance values
mata: __trop_ptd_uval = `ptd_unit_val'
mata: __trop_ptd_tval = `ptd_time_val'

*--- Standardize outcome: (Y - overall mean) / overall SD, in place. ---*
local Ysd = 1
mata: trop_standardize_outcome(Y)

tempvar ever_treated unit_row unit_tag

qui bys `2': egen `ever_treated' = max(`4')
qui egen `unit_row' = group(`2')
qui egen `unit_tag' = tag(`2')

qui levelsof `unit_row' if `ever_treated' == 1 & `unit_tag' == 1, local(treated_rows)
local n_treated : word count `treated_rows'
local ntreated = `n_treated'

qui levelsof `2', local(__trop_uvals)
qui levelsof `3', local(__trop_tvals)

drop `ever_treated' `unit_row' `unit_tag'

if (length("`if'")+length("`in'")>0) {
    restore
}

*------------------------------------------------------------------------------*
* (2) Calculate ATT
*------------------------------------------------------------------------------*
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

* Adaptive search sweeps the grids below literally in phase 1 (so its
* defaults are the joint/cycle grids), then brackets and refines.

if "`unit_grid'" == "" local unit_grid "`default_unit_grid'"
if "`time_grid'" == "" local time_grid "`default_time_grid'"
if "`nn_grid'"   == "" local nn_grid   "`default_nn_grid'"

*K for placebo folds (default 5, overridable via kfold())*
if `kfold' <= 0 local kfold = 5

if `ntrials'  <= 0 local ntrials  = 200

*Run CV if any lambda is unspecified*
local need_cv = (1 - `lambda_unit_set') + (1 - `lambda_time_set') + (1 - `lambda_nn_set')

*--- cv() compatibility checks and ignored-suboption notes ---*
if `need_cv' == 0 & `"`cv'"' != "" {
    di as txt "note: cv() ignored; all three lambdas were fixed."
}
if `need_cv' > 0 {
    if "`group'" == "cell" & inlist("`cv_method'", "resample", "kfold") {
        di as error "cv(`cv_method') is only defined for group(time)."
        di as error "For group(cell), use cv(loocv) or cv(loocv, cells(#))."
        exit 198
    }
    if "`group'" == "time" & "`cv_method'" == "loocv" {
        di as error "cv(loocv) is only defined for group(cell)."
        di as error "For group(time), use cv(resample) or cv(kfold)."
        exit 198
    }
    local _nctrl = `N' - `n_treated'
    if "`cv_method'" == "resample" & "`group'" == "time" & `ntreated' > `_nctrl' {
        di as error "cv(resample): the design has `ntreated' treated units but only `_nctrl' never-treated units."
        di as error "Resample CV stamps one placebo per treated unit, so it cannot run on this design."
        di as error "Use cv(kfold) or fix all three lambdas."
        exit 459
    }

    if "`_trials'" != "" & "`cv_method'" != "resample" ///
        di as txt "note: trials() ignored under cv(`cv_method')."

    if "`_trials'" != "" & "`cv_method'" == "resample" & "`group'" != "time" ///
        di as txt "note: trials() ignored; cv(resample) is only used with group(time)."


    if "`_folds'" != "" & "`cv_method'" != "kfold" ///
        di as txt "note: folds() ignored under cv(`cv_method')."

    if "`_folds'" != "" & "`cv_method'" == "kfold" & "`group'" != "time" ///
        di as txt "note: folds() ignored; cv(kfold) is only used with group(time)."

    if "`_cells'" != "" & !("`group'" == "cell" & "`cv_method'" == "loocv") ///
        di as txt "note: cells() ignored (only used with group(cell) cv(loocv))."

    if "`_points'" != "" & "`cv_search'" != "adaptive" ///
        di as txt "note: points() ignored (only used with cv(..., adaptive))."
    if "`_expansions'" != "" & "`cv_search'" != "adaptive" ///
        di as txt "note: expansions() ignored (only used with cv(..., adaptive))."

    if inlist("`cv_search'", "joint", "adaptive") {
        if `lambda_unit_set' & `"`_unit_grid'"' != "" di as txt "note: unit_grid() ignored; lambda_unit is fixed."
        if `lambda_time_set' & `"`_time_grid'"' != "" di as txt "note: time_grid() ignored; lambda_time is fixed."
        if `lambda_nn_set'   & `"`_nn_grid'"'   != "" di as txt "note: nn_grid() ignored; lambda_nn is fixed."
    }
}

if `need_cv' > 0 {
    if inlist("`cv_search'", "joint", "adaptive") {
        if `lambda_unit_set' local unit_grid "`lu'"
        if `lambda_time_set' local time_grid "`lt'"
        if `lambda_nn_set' {
            if "`lnn'" == "." local nn_grid "."
            else              local nn_grid "`lnn'"
        }
    }

    if "`group'" == "cell" {
        local _ncc = `N' * `T' - `ntrcells'
        if `ncells' > 0 & `ncells' < `_ncc' local _ncclab "`ncells' samples of `_ncc' total"
        else                                local _ncclab "`_ncc'"
        di as txt "Cross-validating lambdas using `cv_method' with `cv_search' search, and `_ncclab' control cells."
        di as txt "To reduce `cv_method' computational time, reduce number of cells or set lambdas."
    }
    else {
        if "`cv_method'" == "resample" {
            di as txt "Cross-validating lambdas using resample with `cv_search' search, and `ntrials' trials (seed `cv_seed')."
            di as txt "To reduce resample computational time, reduce no of trials or set lambdas."
        }
        else {
            di as txt "Cross-validating lambdas using k-fold with `cv_search' search, and `kfold' folds (seed `cv_seed')."
            di as txt "To reduce k-fold computational time, reduce no of folds or set lambdas."
        }
    }

    set seed `cv_seed'

    mata: trop_cv_setup(Y, W)

    if "`verbose'" != "" {
        di as txt ""
        di as txt "Selected lambda: unit, time and nn = [`lu' ; `lt' ; `lnn']"
    }
}

if "`lnn'" == "." local lnn_arg = .
else              local lnn_arg = `lnn'

mata: trop_point_att(Y, W, "`group'", `lu', `lt', `lnn_arg')

* Map point estimates back to raw outcome units
local tau_hat = `tau_hat' * `Ysd'
matrix __trop_ttau = `Ysd' * __trop_ttau

*--- Label heterogeneous groups ---*
local __trop_keys ""
if "`group'" == "cell" {
    forvalues g = 1/`n_targets' {
        local __uv : word `=__trop_tinfo[`g',1]' of `__trop_uvals'
        local __tv : word `=__trop_tinfo[`g',2]' of `__trop_tvals'
        local __trop_keys "`__trop_keys' u`__uv'_t`__tv'"
        matrix __trop_tinfo[`g',1] = `__uv'
        matrix __trop_tinfo[`g',2] = `__tv'
    }
    matrix colnames __trop_tinfo = unit time
    mata: trop_target_grid()
}
else {
    forvalues g = 1/`n_targets' {
        local __sv : word `=__trop_tinfo[`g',1]' of `__trop_tvals'
        local __ev : word `=__trop_tinfo[`g',2]' of `__trop_tvals'
        local __trop_keys "`__trop_keys' t`__sv'_`__ev'"
        matrix __trop_tinfo[`g',1] = `__sv'
        matrix __trop_tinfo[`g',2] = `__ev'
        local __gu ""
        forvalues k = 1/`=colsof(__trop_tunits)' {
            if !missing(__trop_tunits[`g', `k']) {
                local __uv : word `=__trop_tunits[`g',`k']' of `__trop_uvals'
                local __gu "`__gu' `__uv'"
            }
        }
        local __gu = strtrim("`__gu'")
        if `g' == 1 local __trop_units "`__gu'"
        else        local __trop_units "`__trop_units' | `__gu'"
    }
    matrix colnames __trop_tinfo = start end n_units n_cells
    cap matrix drop __trop_tunits
}
cap matrix colnames __trop_ttau  = `__trop_keys'
cap matrix colnames __trop_twt   = `__trop_keys'
cap matrix rownames __trop_tinfo = `__trop_keys'

if `n_zov' > 0 {
    di as txt "warning: `n_zov' comparison unit-group pair(s) had no eligible periods to compute a"
    di as txt "         unit distance from; those units received unit weight 0 in the affected group(s):"
    forvalues k = 1/`n_zov' {
        local __zg = __trop_zovinfo[`k',1]
        local __zr = __trop_zovinfo[`k',2]
        local __zk : word `__zg' of `__trop_keys'
        local __zu : word `__zr' of `__trop_uvals'
        di as txt "         unit `__zu' in group `__zk'"
    }
}
cap matrix drop __trop_zovinfo

*--------------------------------------------------------------------------*
* (3) Standard error: bootstrap or jackknife
*--------------------------------------------------------------------------*
local se = .
local ci_lo = .
local ci_hi = .
local n_jk = .
if "`vce'" == "bootstrap" {
    if `seed' != 0  set seed `seed'
    di as txt "To reduce computational time, reduce reps() or use vce(jackknife)."

    mata: st_numscalar("r(se)", trop_bootstrap_att(Y, W, "`group'", `lu', `lt', `lnn_arg', `reps', `level'))
    local se    = r(se)
    local ci_lo = r(ci_lower)
    local ci_hi = r(ci_upper)
    local se    = `se' * `Ysd'
    local ci_lo = `ci_lo' * `Ysd'
    local ci_hi = `ci_hi' * `Ysd'
}
else if "`vce'" == "jackknife" {
    mata: st_numscalar("r(se)", trop_jackknife_att(Y, W, "`group'", `lu', `lt', `lnn_arg'))
    local se   = r(se)
    local n_jk = r(df_jk)
    local se   = `se' * `Ysd'
    if "`se'" != "." {
        local ci_lo = `tau_hat' - invnormal(1 - (100 - `level')/200) * `se'
        local ci_hi = `tau_hat' + invnormal(1 - (100 - `level')/200) * `se'
    }
}

*--------------------------------------------------------------------------*
* (4) Return output
*--------------------------------------------------------------------------*
ereturn clear
ereturn scalar tau = `tau_hat'
ereturn scalar N = `N'
ereturn scalar T = `T'
ereturn scalar N_treated = `n_treated'
ereturn scalar n_groups = `n_targets'
ereturn local  group "`group'"
ereturn matrix group_tau    = __trop_ttau
ereturn matrix group_weight = __trop_twt
ereturn matrix group_info   = __trop_tinfo
if "`group'" == "time" ereturn local group_units "`__trop_units'"
if "`group'" == "cell" ereturn matrix group_grid = __trop_tgrid
ereturn scalar lambda_unit = `lu'
ereturn scalar lambda_time = `lt'
if "`lnn'" == "." ereturn scalar lambda_nn = .
else              ereturn scalar lambda_nn = `lnn'

if "`group'" == "time" {
    local _ptdu = cond("`ptd_unit_val'" == ".", "default", "`ptd_unit_val'")
    local _ptdt = cond("`ptd_time_val'" == ".", "default", "`ptd_time_val'")
    ereturn local pooled_treat_distance "unit(`_ptdu') time(`_ptdt')"
}

ereturn local cmd "trop"
ereturn local depvar "`1'"
ereturn local unitvar "`unitname'"
ereturn local timevar "`3'"
ereturn local treatvar "`4'"
ereturn local vce "`vce'"

local have_se = ("`se'" != "." & "`se'" != "")
if `have_se' {
    ereturn scalar se = `se'
    ereturn scalar ci_lower = `ci_lo'
    ereturn scalar ci_upper = `ci_hi'
    ereturn scalar z = `tau_hat' / `se'
    ereturn scalar pvalue = 2 * normal(-abs(`tau_hat' / `se'))
    ereturn scalar level = `level'
    if "`vce'" == "jackknife" & "`n_jk'" != "." ereturn scalar N_jack = `n_jk'
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
    di as text %12s " " " {c |}  " as text "(selected by `cv_method' CV)"
}
di as text "{hline 13}{c BT}{hline 50}"

if "`ptd_unit_val'" != "." | "`ptd_time_val'" != "." {
    di as txt "(pooled treated-block distances: unit `=cond("`ptd_unit_val'"==".","default","`ptd_unit_val'")', time `=cond("`ptd_time_val'"==".","default","`ptd_time_val'")')"
}

if `n_targets' > 1 {
    if "`group'" == "cell" di as txt ///
        `"(`n_targets' per-cell effects: {stata "matrix list e(group_grid), format(%9.4f)":e(group_grid)} [unit x time], e(group_tau), e(group_info))"'
    else di as txt ///
        "(`n_targets' per-spell effects in e(group_tau), e(group_weight), e(group_info), e(group_units))"
}

if "`detail'" != "" & `n_targets' > 1 & "`group'" == "time" {
    di as txt ""
    di as txt %12s "time periods" %7s "start" %6s "end" %7s "units" %7s "cells" %11s "tau" "   cohort's units"
    tempname __ti __tt
    matrix `__ti' = e(group_info)
    matrix `__tt' = e(group_tau)
    local __rest `"`e(group_units)'"'
    forvalues g = 1/`n_targets' {
        gettoken __gu __rest : __rest, parse("|")
        if `"`__gu'"' == "|" gettoken __gu __rest : __rest, parse("|")
        local __nm : word `g' of `: colnames `__tt''
        di as txt %12s "`__nm'" as res %7.0g `__ti'[`g',1] %6.0g `__ti'[`g',2] ///
           %7.0g `__ti'[`g',3] %7.0g `__ti'[`g',4] %11.4f `__tt'[1,`g'] ///
           as txt "   `=strtrim(`"`__gu'"')'"
    }
}
end

*------------------------------------------------------------------------------*
* Mata functions
*------------------------------------------------------------------------------*

cap mata: mata drop trop_standardize_outcome()
cap mata: mata drop trop_svt()
cap mata: mata drop trop_cv_single()
cap mata: mata drop trop_cv_cycle()
cap mata: mata drop trop_cv_setup()
cap mata: mata drop trop_kfold_sets()
cap mata: mata drop trop_resample_sets()
cap mata: mata drop trop_cv_joint()
cap mata: mata drop trop_suff_gram()
cap mata: mata drop trop_suff_ginv()
cap mata: mata drop trop_suff_rhs()
cap mata: mata drop trop_coef_to_fit()
cap mata: mata drop trop_nuclear_core()
cap mata: mata drop trop_nuclear_path_suff()
cap mata: mata drop trop_placebo_rmse_path()
cap mata: mata drop trop_time_weights2()
cap mata: mata drop trop_unit_weights2()
cap mata: mata drop trop_target_tau()
cap mata: mata drop trop_build_groups()
cap mata: mata drop trop_att()
cap mata: mata drop trop_point_att()
cap mata: mata drop trop_quantile()
cap mata: mata drop trop_bootstrap_att()
cap mata: mata drop trop_jackknife_att()
cap mata: mata drop trop_ms_now()
cap mata: mata drop trop_ms_fmt()
cap mata: mata drop trop_ms_init()
cap mata: mata drop trop_ms_tick()
cap mata: mata drop trop_loocv_cells()
cap mata: mata drop trop_loocv_cell_rmse_path()
cap mata: mata drop trop_cv_single_cell()
cap mata: mata drop trop_cv_loocv_cell()
cap mata: mata drop trop_cv_joint_cell()
cap mata: mata drop trop_target_grid()
cap mata: mata drop trop_adkey()
cap mata: mata drop trop_adaptive_eval()
cap mata: mata drop trop_cv_adaptive()
cap mata: mata drop trop_ad_rangeline()
cap mata: mata drop trop_ad_heat()
cap mata: mata drop trop_ad_axlist()
cap mata: mata drop trop_ad_num()
cap mata: mata drop trop_ad_gridline()
cap mata: mata drop trop_ad_bracket()
cap mata: mata drop trop_ad_pad()



mata:
void trop_standardize_outcome(real matrix Y)
{
    real scalar m, s, n
    n = rows(Y) * cols(Y)
    m = sum(Y) / n
    s = sqrt( sum((Y :- m):^2) / (n - 1) )
    if (s <= 0 | s >= .) s = 1                 // degenerate / missing -> no rescale
    Y = (Y :- m) :/ s
    st_local("Ysd", strofreal(s, "%21.16g"))
}
end

mata:
real scalar trop_ms_now()
{
    return( clock(c("current_date") + " " + c("current_time"), "DMY hms") )
}
end

mata:
string scalar trop_ms_fmt(real scalar secs)
{
    real scalar s, m, h
    s = floor(secs)
    if (s < 0 | s >= .) s = 0
    h = floor(s / 3600)
    m = floor((s - 3600 * h) / 60)
    s = s - 3600 * h - 60 * m
    if (h > 0) return(sprintf("%g:%02.0f:%02.0f", h, m, s))
    return(sprintf("%g:%02.0f", m, s))
}
end

mata:
// State: (start ms, total units, failures, printed-anything flag, last-line ms).
real rowvector trop_ms_init(real scalar total)
{
    real scalar now
    now = trop_ms_now()
    return( (now, total, 0, 0, now) )
}
end

mata:
void trop_ms_tick(real rowvector S, real scalar i, real scalar ok)
{
    real scalar total, now, el, rem, cur, prev, pct

    total = S[2]
    if (!ok) S[3] = S[3] + 1
    now = trop_ms_now()
    el  = (now - S[1]) / 1000

    if (i == 1 & i < total) {
        if (el >= 2 & el * total >= 10) {
            printf("{txt}  First replication %s, roughly %s expected in total\n",
                   trop_ms_fmt(el), trop_ms_fmt(el * total))
            displayflush()
            S[4] = 1
            S[5] = now
        }
        return
    }
    if (i == total) {
        if (S[4] | el >= 2) {
            printf("{txt}  100%%  (%g/%g)   %s elapsed", total, total, trop_ms_fmt(el))
            if (S[3] > 0) printf("{txt}, %g failed", S[3])
            printf("{txt}\n")
            displayflush()
        }
        return
    }
    cur  = floor(20 * i / total)
    prev = floor(20 * (i - 1) / total)
    if (cur > prev & (now - S[5]) / 1000 >= 30) {
        pct = cur * 5
        rem = el / i * (total - i)
        printf("{txt}  %3.0f%%  (%g/%g)   %s elapsed   about %s remaining\n",
               pct, i, total, trop_ms_fmt(el), trop_ms_fmt(rem))
        displayflush()
        S[4] = 1
        S[5] = now
    }
}
end

mata:
real matrix trop_suff_gram(real matrix delta, real matrix W)
{
    real scalar    N, T, p, S, d, ridge, i
    real colvector ru, cu
    real rowvector ct, cw
    real matrix    A, dW

    N = rows(delta); T = cols(delta); p = N + T
    ridge = 1e-10

    A  = J(p, p, 0)
    S  = sum(delta)
    ru = rowsum(delta)        // N x 1  unit weight sums
    ct = colsum(delta)        // 1 x T  time weight sums

    // intercept row/col (index 1)
    A[1,1] = S
    A[|1,2     \ 1,N      |] = ru[|2\N|]'
    A[|2,1     \ N,1      |] = ru[|2\N|]
    A[|1,N+1   \ 1,N+T-1  |] = ct[|2\T|]
    A[|N+1,1   \ N+T-1,1  |] = ct[|2\T|]'

    // unit-unit and time-time diagonals
    A[|2,2     \ N,N        |] = diag(ru[|2\N|])
    A[|N+1,N+1 \ N+T-1,N+T-1|] = diag(ct[|2\T|])

    // unit-time cross block = delta[unit, time]
    A[|2,N+1   \ N,N+T-1  |] = delta[|2,2 \ N,T|]
    A[|N+1,2   \ N+T-1,N  |] = delta[|2,2 \ N,T|]'

    // treatment column (W is 0/1 so W:^2 = W): weighted sums over treated cells
    dW = delta :* W
    d  = sum(dW)
    cu = rowsum(dW)           // N x 1
    cw = colsum(dW)           // 1 x T
    A[1,p] = d ; A[p,1] = d
    A[|2,p     \ N,p      |] = cu[|2\N|]
    A[|p,2     \ p,N      |] = cu[|2\N|]'
    A[|N+1,p   \ N+T-1,p  |] = cw[|2\T|]'
    A[|p,N+1   \ p,N+T-1  |] = cw[|2\T|]
    A[p,p] = d

    // tiny ridge on the diagonal for numerical safety
    for (i = 1; i <= p; i++) A[i,i] = A[i,i] + ridge

    return(A)
}
end

// Invert the (symmetric, positive-definite) Gram once; reused across
// nuclear-norm path. cholinv() returns all-missing on a near-singular matrix
// so we guard and fall back to the generalized symmetric inverse invsym().
mata:
real matrix trop_suff_ginv(real matrix delta, real matrix W)
{
    real matrix Gi
    Gi = cholinv(trop_suff_gram(delta, W))
    if (missing(Gi)) Gi = invsym(trop_suff_gram(delta, W))
    return(Gi)
}
end

mata:
real colvector trop_suff_rhs(real matrix delta, real matrix Ytil, real matrix W)
{
    real scalar    N, T, p
    real matrix    DYt
    real colvector b, rb
    real rowvector cb

    N = rows(delta); T = cols(delta); p = N + T
    DYt = delta :* Ytil
    rb  = rowsum(DYt)         // N x 1
    cb  = colsum(DYt)         // 1 x T

    b = J(p, 1, 0)
    b[1]            = sum(DYt)
    b[|2\N|]        = rb[|2\N|]
    b[|N+1\N+T-1|]  = cb[|2\T|]'
    b[p]            = sum(DYt :* W)
    return(b)
}
end

mata:
real matrix trop_coef_to_fit(real colvector coef, real matrix W, real scalar N, real scalar T)
{
    real scalar    mu, tau
    real colvector alpha
    real rowvector beta

    mu    = coef[1]
    alpha = 0 \ coef[|2\N|]                 
    beta  = 0 , (coef[|N+1\N+T-1|])'        
    tau   = coef[N+T]
    return( mu :+ (alpha * J(1,T,1)) :+ (J(N,1,1) * beta) :+ tau :* W )
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
void trop_nuclear_core(
    real matrix Y, real matrix W, real matrix delta,
    real scalar lambda_nn, real scalar tol, real scalar max_iter, real scalar cv_mode,
    real matrix Ginv,
    real matrix L, real scalar tau, real scalar iters)
{
    real scalar N, T, p, k, maxd, step, thr
    real scalar a, a_next, mom, nucnorm
    real scalar tau_new, tau_old, obj_new, obj_old, dL, dtau, dobj, normL, restart
    real matrix L_prev, Yk, Ytil, R, Gstep, L_new
    real colvector coef

    N = rows(Y); T = cols(Y); p = N + T
    maxd = max(delta)
    step = 1 / (2 * maxd)
    thr  = step * lambda_nn

    L_prev  = L                 
    a       = 1
    tau_old = .
    obj_old = .
    nucnorm = 0
    iters   = 0

    for (k = 1; k <= max_iter; k++) {
        iters  = k
        a_next = (1 + sqrt(1 + 4 * a^2)) / 2
        mom    = (a - 1) / a_next
        Yk     = L + mom :* (L - L_prev)
        Ytil   = Y - Yk
        coef   = Ginv * trop_suff_rhs(delta, Ytil, W)
        R      = Ytil - trop_coef_to_fit(coef, W, N, T)
        Gstep  = Yk + (delta :* R) :/ maxd
        L_new  = trop_svt(Gstep, thr, nucnorm)

        restart = sum((Yk - L_new) :* (L_new - L))
        if (restart > 0) a_next = 1

        if (cv_mode) {
            tau_new = coef[p]
            dtau = (tau_old < . ? abs(tau_new - tau_old) / (1 + abs(tau_old)) : .)
            L_prev = L; L = L_new; a = a_next; tau_old = tau_new
            if (dtau < tol) break
        }
        else {
            coef    = Ginv * trop_suff_rhs(delta, Y - L_new, W)
            R       = (Y - L_new) - trop_coef_to_fit(coef, W, N, T)
            obj_new = sum(delta :* (R:^2)) + lambda_nn * nucnorm
            tau_new = coef[p]
            normL   = sqrt(sum(L:^2))
            dL   = sqrt(sum((L_new - L):^2)) / (1 + normL)
            dtau = (tau_old < . ? abs(tau_new - tau_old) / (1 + abs(tau_old)) : .)
            dobj = (obj_old < . ? abs(obj_new - obj_old) / (1 + abs(obj_old)) : .)
            L_prev = L; L = L_new; a = a_next; tau_old = tau_new; obj_old = obj_new
            if (dtau < tol & dobj < tol) break
        }
    }
    coef = Ginv * trop_suff_rhs(delta, Y - L, W)
    tau  = coef[p]
}
end

mata:
real rowvector trop_nuclear_path_suff(
    real matrix Y, real matrix W, real matrix delta,
    real rowvector nn_grid, real scalar tol, real scalar max_iter, real scalar cv_mode,
    real rowvector iters_out)
{
    real scalar ng, p, nfin, k, j, pos, lam, tau, it, twls, i
    real matrix Ginv, Lw
    real rowvector taus, finvals, fin_pos, inf_pos
    real colvector ord, coef

    ng = cols(nn_grid)
    p  = rows(Y) + cols(Y)
    taus      = J(1, ng, .)
    iters_out = J(1, ng, 0)

    Ginv = trop_suff_ginv(delta, W)         // built ONCE, reused across the path

    coef = Ginv * trop_suff_rhs(delta, Y, W)   // WLS tau for any inf (.) entries
    twls = coef[p]

    fin_pos = selectindex(nn_grid :< .)     // positions of finite lambda_nn (row in -> row out)
    inf_pos = selectindex(nn_grid :>= .)    // positions of infinite (missing) lambda_nn

    if (cols(fin_pos) > 0) {
        nfin    = cols(fin_pos)
        finvals = nn_grid[fin_pos]
        ord     = order(finvals', 1)        // ascending; walked descending below
        Lw      = J(rows(Y), cols(Y), 0)    // cold start for the LARGEST lambda_nn
        for (k = nfin; k >= 1; k--) {
            j   = ord[k]
            lam = finvals[j]
            pos = fin_pos[j]
            tau = .; it = .
            trop_nuclear_core(Y, W, delta, lam, tol, max_iter, cv_mode, Ginv, Lw, tau, it)
            taus[pos]      = tau
            iters_out[pos] = it             // Lw now warm-starts the next smaller lambda_nn
        }
    }
    for (i = 1; i <= cols(inf_pos); i++) {
        taus[inf_pos[i]]      = twls
        iters_out[inf_pos[i]] = 0
    }
    return(taus)
}
end

mata:
real rowvector trop_placebo_rmse_path(
    real matrix Yc,
    pointer(real colvector) rowvector sets,
    real scalar lu, real scalar lt, real rowvector nn_grid,
    real scalar tol, real scalar max_iter, real scalar cv_mode)
{
    real scalar Nc, T, J_sets, j, ng, g, G, q, np, sz, wsum, g2
    real colvector set_rows, us, omega
    real rowvector per, theta, taus_g, it_out
    real rowvector sq_accum, count, rmse, att_curve, badnn
    real matrix Wp, gid, mask, delta, f
    external real matrix __trop_Pat

    Nc = rows(Yc); T = cols(Yc); J_sets = cols(sets); ng = cols(nn_grid)
    np = rows(__trop_Pat)
    if (np < 1) {
        errprintf("trop_placebo_rmse_path(): empty pattern bank.\n")
        exit(459)
    }
    sq_accum = J(1, ng, 0)
    count    = J(1, ng, 0)

    for (j = 1; j <= J_sets; j++) {
        set_rows = *sets[j]
        Wp = J(Nc, T, 0)
        for (q = 1; q <= rows(set_rows); q++)
            Wp[set_rows[q], .] = __trop_Pat[mod(q - 1, np) + 1, .]

        G = .
        gid = trop_build_groups(Wp, "time", G)
        att_curve = J(1, ng, 0)
        badnn     = J(1, ng, 0)
        wsum = 0
        for (g = 1; g <= G; g++) {
            mask  = (gid :== g)
            us    = selectindex(rowsum(mask) :> 0)
            per   = selectindex(colsum(mask) :> 0)
            sz    = sum(mask)
            omega = trop_unit_weights2(Yc, Wp, us, per, lu, 1)
            theta = trop_time_weights2(T, per, lt, 1)
            f     = (1 :- Wp) :+ Wp :* mask
            delta = f :* (omega * theta)
            it_out = .
            taus_g = trop_nuclear_path_suff(Yc, mask, delta, nn_grid,
                                            tol, max_iter, cv_mode, it_out)
            for (g2 = 1; g2 <= ng; g2++) {
                if (taus_g[g2] < .) att_curve[g2] = att_curve[g2] + sz * taus_g[g2]
                else                badnn[g2]     = badnn[g2] + 1
            }
            wsum = wsum + sz
        }
        att_curve = att_curve :/ wsum
        for (g = 1; g <= ng; g++) {
            if (badnn[g] == 0) {
                sq_accum[g] = sq_accum[g] + att_curve[g]^2
                count[g]    = count[g] + 1
            }
        }
    }

    rmse = J(1, ng, .)
    for (g = 1; g <= ng; g++) {
        if (count[g] > 0) rmse[g] = sqrt(sq_accum[g] / count[g])
    }
    return(rmse)
}
end

mata:
real scalar trop_cv_single(
    real matrix Yc,
    pointer(real colvector) rowvector sets,
    real rowvector grid,
    real scalar which_lambda,
    real scalar fixed1,
    real scalar fixed2,
    real scalar tol,
    real scalar max_iter,
    real rowvector scores,
    string scalar lab,
    real scalar show
)
{
    real scalar ng, g, best_idx, best_score, cv_mode
    real rowvector nn1

    ng = cols(grid)
    cv_mode = 1                              // CV uses the fast tau-stop solver

    if (which_lambda == 3) {
        // nn coordinate: lu=fixed1, lt=fixed2; whole grid in ONE warm sweep/set
        scores = trop_placebo_rmse_path(Yc, sets,
                                        fixed1, fixed2, grid, tol, max_iter, cv_mode)
    }
    else if (which_lambda == 1) {
        // unit coordinate: lu=grid[g], lt=fixed1, lnn=fixed2
        scores = J(1, ng, .)
        for (g = 1; g <= ng; g++) {
            nn1 = (fixed2)
            scores[g] = trop_placebo_rmse_path(Yc, sets,
                                               grid[g], fixed1, nn1, tol, max_iter, cv_mode)[1]
        }
    }
    else if (which_lambda == 2) {
        // time coordinate: lu=fixed1, lt=grid[g], lnn=fixed2
        scores = J(1, ng, .)
        for (g = 1; g <= ng; g++) {
            nn1 = (fixed2)
            scores[g] = trop_placebo_rmse_path(Yc, sets,
                                               fixed1, grid[g], nn1, tol, max_iter, cv_mode)[1]
        }
    }
    else {
        errprintf("trop_cv_single(): which_lambda must be 1, 2, or 3.\n")
        exit(198)
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
    if (show) {
        printf("{txt}%s -> %s   ", lab,
               (grid[best_idx] >= . ? "inf" : sprintf("%g", grid[best_idx])))
        displayflush()
    }
    return(grid[best_idx])
}
end

mata:
real scalar trop_cv_cycle(
    real matrix Yc,
    pointer(real colvector) rowvector sets,
    real rowvector unit_grid, real rowvector time_grid, real rowvector nn_grid,
    real scalar tune_unit, real scalar tune_time, real scalar tune_nn,
    real scalar tol, real scalar max_iter, real scalar max_cycles,
    real scalar lambda_unit, real scalar lambda_time, real scalar lambda_nn,
    real scalar cycles_used)
{
    real scalar c, lu_old, lt_old, lnn_old, t0
    real rowvector scores, fin

    t0 = trop_ms_now()

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

        printf("{txt}  cycle %g of up to %g:  ", c, max_cycles)
        displayflush()
        if (tune_unit)
            lambda_unit = trop_cv_single(Yc, sets, unit_grid, 1,
                                         lambda_time, lambda_nn, tol, max_iter, scores,
                                         "lambda_unit", 1)
        if (tune_time)
            lambda_time = trop_cv_single(Yc, sets, time_grid, 2,
                                         lambda_unit, lambda_nn, tol, max_iter, scores,
                                         "lambda_time", 1)
        if (tune_nn)
            lambda_nn = trop_cv_single(Yc, sets, nn_grid, 3,
                                       lambda_unit, lambda_time, tol, max_iter, scores,
                                       "lambda_nn", 1)
        printf("{txt}%s elapsed\n", trop_ms_fmt((trop_ms_now() - t0) / 1000))
        displayflush()

        if (lambda_unit == lu_old & lambda_time == lt_old & lambda_nn == lnn_old) {
            cycles_used = c
            printf("{txt}  converged after %g cycle(s), %s total\n",
                   c, trop_ms_fmt((trop_ms_now() - t0) / 1000))
            displayflush()
            return(lambda_unit)
        }
    }
    cycles_used = max_cycles
    printf("{txt}  reached max cycles (%g) without full convergence, %s total\n",
           max_cycles, trop_ms_fmt((trop_ms_now() - t0) / 1000))
    displayflush()
    return(lambda_unit)
}
end

mata:
void trop_cv_setup(real matrix Yfull, real matrix Wfull)
{
    real matrix Yc, cells
    real rowvector unit_grid, time_grid, nn_grid
    real colvector control_rows
    external real matrix __trop_Pat
    real scalar K, tol, max_iter, max_cycles, Nc
    real scalar tune_unit, tune_time, tune_nn
    real scalar lambda_unit, lambda_time, lambda_nn, cycles_used
    real scalar ntrials, ntreated, best_rmse, n_eval
    real scalar nn_inf_only, npts_a, nexp_a, phases_a
    real rowvector fin_r, res_a
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

    __trop_Pat = Wfull[selectindex(rowsum(Wfull) :> 0), .]

    unit_grid = strtoreal(tokens(st_local("unit_grid")))
    time_grid = strtoreal(tokens(st_local("time_grid")))
    nn_grid   = strtoreal(tokens(st_local("nn_grid")))

    K          = strtoreal(st_local("kfold"))
    tol        = 1e-6
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

    cycles_used = .

    // Adaptive search: phase 1 sweeps the grids above literally (the defaults
    // are the joint/cycle grids); later phases bracket the winner between its
    // neighbouring grid points and refine on uniform points() grids.
    if (search == "adaptive") {
        fin_r = select(nn_grid, nn_grid :< .)
        nn_inf_only = (cols(fin_r) == 0)        // fixed at inf (or all-. grid)
        npts_a = strtoreal(st_local("npoints"))
        nexp_a = strtoreal(st_local("nexpansions"))
    }

    if (st_local("group") == "cell") {
    cells = trop_loocv_cells(Wfull, strtoreal(st_local("ncells")))
    if (search == "adaptive") {
            phases_a = .
            res_a = trop_cv_adaptive(Yfull, Wfull, J(1, 0, NULL), cells, 1,
                                     unit_grid, tune_unit,
                                     time_grid, tune_time,
                                     nn_grid, tune_nn, nn_inf_only,
                                     npts_a, nexp_a, tol, max_iter,
                                     (st_local("verbose") != ""), phases_a)
            lambda_unit = res_a[1]; lambda_time = res_a[2]; lambda_nn = res_a[3]
            st_local("cv_cycles", strofreal(phases_a))
        }
        else if (search == "joint") {
            best_rmse = .; n_eval = .
            trop_cv_joint_cell(Yfull, Wfull, cells,
                               unit_grid, time_grid, nn_grid,
                               tol, max_iter,
                               lambda_unit, lambda_time, lambda_nn,
                               best_rmse, n_eval)
            st_local("cv_cycles", "0")
        }
        else {
            trop_cv_loocv_cell(Yfull, Wfull, cells,
                               unit_grid, time_grid, nn_grid,
                               tune_unit, tune_time, tune_nn,
                               tol, max_iter, max_cycles,
                               lambda_unit, lambda_time, lambda_nn, cycles_used)
            st_local("cv_cycles", strofreal(cycles_used))
        }
    }
    else {
        // Build placebo sets ONCE, by the chosen sampling method. The command
        // body has already `set seed`'d, so these draws are reproducible.
        if (method == "resample") sets = trop_resample_sets(Nc, ntrials, ntreated)
        else                      sets = trop_kfold_sets(Nc, K)

        if (search == "adaptive") {
            phases_a = .
            res_a = trop_cv_adaptive(Yc, J(0, 0, .), sets, J(0, 2, .), 0,
                                     unit_grid, tune_unit,
                                     time_grid, tune_time,
                                     nn_grid, tune_nn, nn_inf_only,
                                     npts_a, nexp_a, tol, max_iter,
                                     (st_local("verbose") != ""), phases_a)
            lambda_unit = res_a[1]; lambda_time = res_a[2]; lambda_nn = res_a[3]
            st_local("cv_cycles", strofreal(phases_a))
        }
        else if (search == "joint") {
            best_rmse = .; n_eval = .
            trop_cv_joint(Yc, sets,
                          unit_grid, time_grid, nn_grid,
                          tol, max_iter,
                          lambda_unit, lambda_time, lambda_nn, best_rmse, n_eval)
            st_local("cv_cycles", "0")   // joint is not iterative
        }
        else {   // "cycle"
            (void) trop_cv_cycle(Yc, sets,
                                 unit_grid, time_grid, nn_grid,
                                 tune_unit, tune_time, tune_nn,
                                 tol, max_iter, max_cycles,
                                 lambda_unit, lambda_time, lambda_nn, cycles_used)
            st_local("cv_cycles", strofreal(cycles_used))
        }
    }

    st_local("lu", strofreal(lambda_unit))
    st_local("lt", strofreal(lambda_time))
    if (lambda_nn >= .) st_local("lnn", ".")
    else                st_local("lnn", strofreal(lambda_nn))
}
end

mata:
pointer(real colvector) rowvector trop_kfold_sets(real scalar Nc, real scalar K)
{
    real colvector perm
    real scalar base, rem, f, start, len
    pointer(real colvector) rowvector sets

    if (K > Nc) K = Nc                      // avoid empty folds / bad subscripts
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
    real matrix Yc,
    pointer(real colvector) rowvector sets,
    real rowvector unit_grid, real rowvector time_grid, real rowvector nn_grid,
    real scalar tol, real scalar max_iter,
    real scalar lambda_unit, real scalar lambda_time, real scalar lambda_nn,
    real scalar best_rmse, real scalar n_eval)
{
    real scalar    nu, nt, ng, ii, jj, g, t0
    real scalar    best_u, best_t, best_n
    real rowvector scores, S

    nu = cols(unit_grid); nt = cols(time_grid); ng = cols(nn_grid)
    best_rmse = .; best_u = .; best_t = .; best_n = .
    n_eval = 0

    t0 = trop_ms_now()
    printf("{txt}  joint search: evaluating %g (unit x time) pairs, %g nn values each\n",
           nu * nt, ng)
    displayflush()
    S = trop_ms_init(nu * nt)
    for (ii = 1; ii <= nu; ii++) {
        for (jj = 1; jj <= nt; jj++) {
            scores = trop_placebo_rmse_path(Yc, sets,
                         unit_grid[ii], time_grid[jj], nn_grid,
                         tol, max_iter, 1)
            n_eval = n_eval + ng
            for (g = 1; g <= ng; g++) {
                if (scores[g] < . & (best_rmse >= . | scores[g] < best_rmse)) {
                    best_rmse = scores[g]
                    best_u    = unit_grid[ii]
                    best_t    = time_grid[jj]
                    best_n    = nn_grid[g]      // may be missing (.) -> inf/WLS, correct
                }
            }
            trop_ms_tick(S, (ii - 1) * nt + jj, 1)
        }
    }

    if (best_rmse >= .) {
        errprintf("trop_cv_joint(): all CV scores are missing.\n")
        exit(498)
    }

    lambda_unit = best_u
    lambda_time = best_t
    lambda_nn   = best_n
    printf("{txt}  best (%s, %s, %s), %s total\n",
           sprintf("%g", best_u), sprintf("%g", best_t),
           (best_n >= . ? "inf" : sprintf("%g", best_n)),
           trop_ms_fmt((trop_ms_now() - t0) / 1000))
    displayflush()
}
end

mata:
real rowvector trop_time_weights2(real scalar T, real rowvector tp1,
                                  real scalar lt, real scalar pooled)
{
    real rowvector s0, tp0, d
    real scalar    c
    external real scalar __trop_ptd_tval
    if (cols(tp1) == 0 | lt == 0) return(J(1, T, 1))
    s0  = (0..(T-1))
    tp0 = tp1 :- 1
    if (pooled) {
        c = (min(tp0) + max(tp0)) / 2
        d = abs(s0 :- c)
        if (__trop_ptd_tval < .) d[tp1] = J(1, cols(tp1), __trop_ptd_tval)   // literal in-spell dist
    }
    else {
        d = abs(s0 :- tp0[1])
    }
    return(exp(-lt :* d))
}
end

mata:
real colvector trop_unit_weights2(real matrix Y, real matrix W,
    real colvector tu1, real rowvector tp1, real scalar lu, real scalar pooled)
{
    real scalar    N, T, i
    real rowvector avg, maskt, Wi
    real colvector den, d, ismem, zov
    real matrix    M, sqdiff
    external real scalar __trop_ptd_uval
    external real colvector __trop_zov

    N = rows(Y); T = cols(Y)
    __trop_zov = J(0, 1, .)
    if (lu == 0) return(J(N, 1, 1))
    ismem = J(N, 1, 0)
    ismem[tu1] = J(rows(tu1), 1, 1)
    if (pooled) {
        maskt  = (colsum(W[tu1, .]) :== 0)              // no group member treated at s
        maskt[1, tp1] = J(1, cols(tp1), 0)              // exclude target block 
        avg    = mean(Y[tu1, .])                     
        sqdiff = (J(N,1,1) * avg :- Y) :^ 2
        if (__trop_ptd_uval < .) {                     
            sqdiff = sqdiff :* (1 :- ismem) :+ ismem :* (__trop_ptd_uval^2)
        }
        M = (J(N,1,1) * maskt) :* (1 :- W)              // exclude any treated unit j 
    }
    else {
        i      = tu1[1]
        Wi     = W[i, .]
        Wi[1, tp1] = J(1, cols(tp1), 1)                 // 1{u != t}: drop target period
        M      = (J(N,1,1) * (1 :- Wi)) :* (1 :- W)   
        sqdiff = (J(N,1,1) * Y[i, .] :- Y) :^ 2
    }
    den = rowsum(M)
    zov = (den :== 0) :* (1 :- ismem)                  
    __trop_zov = selectindex(zov)
    den = den :+ (den :== 0)                          
    d   = sqrt(rowsum(M :* sqdiff) :/ den)
    return( exp(-lu :* d) :* (1 :- zov) )            
}
end


mata:
real scalar trop_target_tau(real matrix Y, real matrix W, real matrix target_mask,
    real colvector tu1, real rowvector tp1,
    real scalar lu, real scalar lt, real scalar lnn,
    real scalar pooled, real scalar cv_mode)
{
    real scalar    T, mxit
    real colvector omega
    real rowvector theta, taus, it_out
    real matrix    f, delta

    T = cols(Y)
    omega = trop_unit_weights2(Y, W, tu1, tp1, lu, pooled)
    theta = trop_time_weights2(T, tp1, lt, pooled)
    f     = (1 :- W) :+ W :* target_mask
    delta = f :* (omega * theta)

    mxit   = (cv_mode ? 3000 : 5000)
    it_out = .
    taus   = trop_nuclear_path_suff(Y, target_mask, delta, (lnn),
                                    1e-6, mxit, cv_mode, it_out)
    return(taus[1])
}
end

mata:
real matrix trop_build_groups(real matrix W, string scalar grp, real scalar G)
{
    real scalar N, T, i, t, a, b, g, nspell, k, found
    real matrix gid, spells, keys

    N = rows(W); T = cols(W)
    gid = J(N, T, 0)

    if (grp == "cell") {
        g = 0
        for (i = 1; i <= N; i++) {
            for (t = 1; t <= T; t++) {
                if (W[i,t] > 0) {
                    g++
                    gid[i,t] = g
                }
            }
        }
        G = g
        return(gid)
    }

    // grp == "time"
    spells = J(0, 3, .)                         // rows: (start, end, unit)
    for (i = 1; i <= N; i++) {
        t = 1
        while (t <= T) {
            if (W[i,t] > 0) {
                a = t
                while (t <= T) {
                    if (W[i,t] == 0) break
                    t++
                }
                b = t - 1
                spells = spells \ (a, b, i)
            }
            else {
                t++
            }
        }
    }
    nspell = rows(spells)
    G = 0
    keys = J(0, 2, .)                           // unique (start,end)
    for (k = 1; k <= nspell; k++) {
        a = spells[k,1]; b = spells[k,2]; i = spells[k,3]
        found = 0
        for (g = 1; g <= G; g++) {
            if (keys[g,1] == a & keys[g,2] == b) {
                found = g
                break
            }
        }
        if (found == 0) {
            G++
            keys = keys \ (a, b)
            found = G
        }
        gid[i, (a..b)] = J(1, b - a + 1, found)
    }
    return(gid)
}
end

mata:
real scalar trop_att(real matrix Y, real matrix W, real matrix gid, real scalar G,
    real scalar lu, real scalar lt, real scalar lnn,
    real scalar pooled, real scalar cv_mode,
    real rowvector taus_out, real rowvector wts_out)
{
    real scalar    N, T, g, k, show
    real matrix    mask
    real colvector us
    real rowvector per, sizes, S
    external real colvector __trop_zov
    external real matrix    __trop_zov_all

    N = rows(Y); T = cols(Y)
    taus_out = J(1, G, .)
    sizes    = J(1, G, .)
    show = (cv_mode == 0 & G >= 50)
    S = J(1, 5, .)
    if (show) {
        printf("{txt}Computing %g group effects.\n", G)
        displayflush()
        S = trop_ms_init(G)
    }
    for (g = 1; g <= G; g++) {
        mask = (gid :== g)
        us   = selectindex(rowsum(mask) :> 0)   
        per  = selectindex(colsum(mask) :> 0)  
        sizes[g]    = sum(mask)
        taus_out[g] = trop_target_tau(Y, W, mask, us, per, lu, lt, lnn, pooled, cv_mode)
        if (cv_mode == 0) {                    
            for (k = 1; k <= rows(__trop_zov); k++)
                __trop_zov_all = __trop_zov_all \ (g, __trop_zov[k])
        }
        if (show) trop_ms_tick(S, g, taus_out[g] < .)
    }
    wts_out = sizes :/ sum(sizes)
    return( sum(wts_out :* taus_out) )
}
end

mata:
void trop_point_att(real matrix Y, real matrix W, string scalar grp,
                    real scalar lu, real scalar lt, real scalar lnn)
{
    real scalar    G, att, pooled, g, k, maxU
    real matrix    gid, info, mask, Umat
    real rowvector taus, wts, per
    real colvector us
    external real matrix __trop_zov_all

    G = .
    gid    = trop_build_groups(W, grp, G)
    pooled = (grp == "time")
    taus = .; wts = .
    __trop_zov_all = J(0, 2, .)
    att = trop_att(Y, W, gid, G, lu, lt, lnn, pooled, 0, taus, wts)

    if (grp == "cell") {
        info = J(G, 2, .)
        for (g = 1; g <= G; g++) {
            mask = (gid :== g)
            info[g,1] = selectindex(rowsum(mask) :> 0)   // unit row index
            info[g,2] = selectindex(colsum(mask) :> 0)   // time index
        }
    }
    else {
        maxU = 0
        for (g = 1; g <= G; g++) {
            us = selectindex(rowsum(gid :== g) :> 0)
            maxU = max((maxU, rows(us)))
        }
        info = J(G, 4, .)
        Umat = J(G, maxU, .)
        for (g = 1; g <= G; g++) {
            mask = (gid :== g)
            us  = selectindex(rowsum(mask) :> 0)
            per = selectindex(colsum(mask) :> 0)
            info[g,1] = per[1]                            // spell start index
            info[g,2] = per[cols(per)]                    // spell end index
            info[g,3] = rows(us)                          // units in spell
            info[g,4] = sum(mask)                         // treated cells
            for (k = 1; k <= rows(us); k++)               // member unit row
                Umat[g, k] = us[k]                        //   indices, padded
        }
        st_matrix("__trop_tunits", Umat)                  // temp transfer only
    }

    st_local("tau_hat",   strofreal(att, "%21.16g"))
    st_local("n_targets", strofreal(G))
    st_local("n_zov",     strofreal(rows(__trop_zov_all)))
    st_matrix("__trop_zovinfo", __trop_zov_all)
    st_matrix("__trop_ttau", taus)
    st_matrix("__trop_twt",  wts)
    st_matrix("__trop_tinfo", info)
}
end

mata:
real scalar trop_quantile(real colvector x, real scalar q)
{
    real colvector xs
    real scalar    n, h, lo
    xs = sort(x, 1)
    n  = rows(xs)
    if (n == 1) return(xs[1])
    h  = (n - 1) * q + 1
    lo = floor(h)
    if (lo >= n) return(xs[n])
    if (lo < 1)  return(xs[1])
    return( xs[lo] + (h - lo) * (xs[lo+1] - xs[lo]) )
}
end

mata:
real scalar trop_bootstrap_att(
    real matrix Y, real matrix W, string scalar grp,
    real scalar lu, real scalar lt, real scalar lnn,
    real scalar reps, real scalar level)
{
    real scalar    b, att_b, G, pooled, B_eff, a, n1, n0
    real colvector trt_rows, ctrl_rows, samp_trt, samp_ctrl, samp_rows, boot
    real rowvector taus, wts, S
    real matrix    Yb, Wb, gid

    trt_rows  = selectindex(rowsum(W) :> 0)
    ctrl_rows = selectindex(rowsum(W) :== 0)
    n1 = rows(trt_rows); n0 = rows(ctrl_rows)
    if (n1 == 0 | n0 == 0) {
        errprintf("Bootstrap requires both treated and control units.\n")
        return(.)
    }
    pooled = (grp == "time")
    boot = J(reps, 1, .)
    printf("{txt}Bootstrap inference using %g bootstrap replications.\n", reps)
    displayflush()
    S = trop_ms_init(reps)
    for (b = 1; b <= reps; b++) {
        samp_trt  = trt_rows[ ceil(n1 :* runiform(n1, 1)) ]
        samp_ctrl = ctrl_rows[ ceil(n0 :* runiform(n0, 1)) ]
        samp_rows = samp_trt \ samp_ctrl
        Yb = Y[samp_rows, .]
        Wb = W[samp_rows, .]
        if (sum(Wb) == 0) {
            trop_ms_tick(S, b, 0)
            continue
        }
        G = .
        gid  = trop_build_groups(Wb, grp, G)
        taus = .; wts = .
        att_b = trop_att(Yb, Wb, gid, G, lu, lt, lnn, pooled, 1, taus, wts)
        boot[b] = att_b
        trop_ms_tick(S, b, boot[b] < .)
    }
    boot  = select(boot, boot :< .)
    B_eff = rows(boot)
    if (B_eff < 2) {
        errprintf("Too few finite bootstrap replications for an SE.\n")
        return(.)
    }
    a = (100 - level) / 2
    st_numscalar("r(ci_lower)", trop_quantile(boot, a/100))
    st_numscalar("r(ci_upper)", trop_quantile(boot, 1 - a/100))
    return( sqrt( sum((boot :- mean(boot)):^2) / (B_eff - 1) ) )
}
end

mata:
real scalar trop_jackknife_att(
    real matrix Y, real matrix W, string scalar grp,
    real scalar lu, real scalar lt, real scalar lnn)
{
    real scalar    i, N, G, pooled, n_eff, m, att_i
    real colvector keep, loo
    real rowvector taus, wts, S
    real matrix    Yi, Wi, gid

    N = rows(Y)
    pooled = (grp == "time")
    loo = J(N, 1, .)
    printf("{txt}Jackknife inference using %g leave-one-out replications.\n", N)
    displayflush()
    S = trop_ms_init(N)
    for (i = 1; i <= N; i++) {
        if (i == 1)      keep = (2::N)
        else if (i == N) keep = (1::N-1)
        else             keep = (1::i-1) \ ((i+1)::N)
        Yi = Y[keep, .]
        Wi = W[keep, .]
        if (sum(Wi) == 0) {                    // deleted the only treated unit
            trop_ms_tick(S, i, 0)
            continue
        }
        G = .
        gid  = trop_build_groups(Wi, grp, G)
        taus = .; wts = .
        att_i = trop_att(Yi, Wi, gid, G, lu, lt, lnn, pooled, 1, taus, wts)
        loo[i] = att_i
        trop_ms_tick(S, i, att_i < .)
    }
    loo   = select(loo, loo :< .)
    n_eff = rows(loo)
    if (n_eff < 2) {
        errprintf("Too few finite jackknife replications for an SE.\n")
        return(.)
    }
    m = mean(loo)
    st_numscalar("r(df_jk)", n_eff)
    return( sqrt( (n_eff - 1) / n_eff * sum((loo :- m):^2) ) )
}
end

mata:
real matrix trop_loocv_cells(real matrix W, real scalar ncells)
{
    real scalar    N, T, i, t, k, n0
    real matrix    cells
    real colvector keep

    N = rows(W); T = cols(W)
    n0 = sum(W :== 0)
    cells = J(n0, 2, .)
    k = 0
    for (i = 1; i <= N; i++) {
        for (t = 1; t <= T; t++) {
            if (W[i,t] == 0) {
                k++
                cells[k, .] = (i, t)
            }
        }
    }
    if (ncells > 0 & ncells < n0) {
        keep  = order(runiform(n0, 1), 1)[|1 \ ncells|]
        cells = cells[keep, .]
    }
    return(cells)
}
end

mata:
real rowvector trop_loocv_cell_rmse_path(
    real matrix Y, real matrix W, real matrix cells,
    real scalar lu, real scalar lt, real rowvector nn_grid,
    real scalar tol, real scalar max_iter)
{
    real scalar    N, T, k, c, g, ng, i, t
    real matrix    mask, delta, f
    real colvector omega
    real rowvector theta, taus, it_out, sq, cnt, rmse

    N = rows(Y); T = cols(Y); k = rows(cells); ng = cols(nn_grid)
    sq  = J(1, ng, 0)
    cnt = J(1, ng, 0)
    for (c = 1; c <= k; c++) {
        i = cells[c,1]; t = cells[c,2]
        mask = J(N, T, 0)
        mask[i,t] = 1
        omega = trop_unit_weights2(Y, W, i, t, lu, 0)
        theta = trop_time_weights2(T, t, lt, 0)
        f     = (1 :- W) :+ W :* mask
        delta = f :* (omega * theta)
        it_out = .
        taus = trop_nuclear_path_suff(Y, mask, delta, nn_grid, tol, max_iter, 1, it_out)
        for (g = 1; g <= ng; g++) {
            if (taus[g] < .) {
                sq[g]  = sq[g] + taus[g]^2
                cnt[g] = cnt[g] + 1
            }
        }
    }
    rmse = J(1, ng, .)
    for (g = 1; g <= ng; g++) {
        if (cnt[g] > 0) rmse[g] = sqrt(sq[g] / cnt[g])
    }
    return(rmse)
}
end

mata:
real scalar trop_cv_single_cell(
    real matrix Y, real matrix W, real matrix cells,
    real rowvector grid, real scalar which_lambda,
    real scalar fixed1, real scalar fixed2,
    real scalar tol, real scalar max_iter, real rowvector scores,
    string scalar lab, real scalar show)
{
    real scalar    ng, g, best_idx, best_score
    real rowvector nn1

    ng = cols(grid)
    if (which_lambda == 3) {
        scores = trop_loocv_cell_rmse_path(Y, W, cells, fixed1, fixed2, grid,
                                           tol, max_iter)
    }
    else if (which_lambda == 1) {
        scores = J(1, ng, .)
        for (g = 1; g <= ng; g++) {
            nn1 = (fixed2)
            scores[g] = trop_loocv_cell_rmse_path(Y, W, cells, grid[g], fixed1,
                                                  nn1, tol, max_iter)[1]
        }
    }
    else if (which_lambda == 2) {
        scores = J(1, ng, .)
        for (g = 1; g <= ng; g++) {
            nn1 = (fixed2)
            scores[g] = trop_loocv_cell_rmse_path(Y, W, cells, fixed1, grid[g],
                                                  nn1, tol, max_iter)[1]
        }
    }
    else {
        errprintf("trop_cv_single_cell(): which_lambda must be 1, 2, or 3.\n")
        exit(198)
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
        errprintf("trop_cv_single_cell(): all CV scores are missing.\n")
        exit(498)
    }
    if (show) {
        printf("{txt}%s -> %s   ", lab,
               (grid[best_idx] >= . ? "inf" : sprintf("%g", grid[best_idx])))
        displayflush()
    }
    return(grid[best_idx])
}
end

mata:
void trop_cv_loocv_cell(
    real matrix Y, real matrix W, real matrix cells,
    real rowvector unit_grid, real rowvector time_grid, real rowvector nn_grid,
    real scalar tune_unit, real scalar tune_time, real scalar tune_nn,
    real scalar tol, real scalar max_iter, real scalar max_cycles,
    real scalar lambda_unit, real scalar lambda_time, real scalar lambda_nn,
    real scalar cycles_used)
{
    real scalar    c, lu_old, lt_old, lnn_old, t0
    real rowvector scores

    t0 = trop_ms_now()
    scores = .
    printf("{txt}  marginal:  ")
    displayflush()
    if (tune_time) lambda_time = trop_cv_single_cell(Y, W, cells, time_grid, 2,
                                     0, ., tol, max_iter, scores, "lambda_time", 1)
    if (tune_unit) lambda_unit = trop_cv_single_cell(Y, W, cells, unit_grid, 1,
                                     0, ., tol, max_iter, scores, "lambda_unit", 1)
    if (tune_nn)   lambda_nn   = trop_cv_single_cell(Y, W, cells, nn_grid, 3,
                                     0, 0, tol, max_iter, scores, "lambda_nn", 1)
    printf("{txt}%s elapsed\n", trop_ms_fmt((trop_ms_now() - t0) / 1000))
    displayflush()

    for (c = 1; c <= max_cycles; c++) {
        lu_old = lambda_unit; lt_old = lambda_time; lnn_old = lambda_nn
        printf("{txt}  cycle %g of up to %g:  ", c, max_cycles)
        displayflush()
        if (tune_unit) lambda_unit = trop_cv_single_cell(Y, W, cells, unit_grid, 1,
                                         lambda_time, lambda_nn, tol, max_iter, scores,
                                         "lambda_unit", 1)
        if (tune_time) lambda_time = trop_cv_single_cell(Y, W, cells, time_grid, 2,
                                         lambda_unit, lambda_nn, tol, max_iter, scores,
                                         "lambda_time", 1)
        if (tune_nn)   lambda_nn   = trop_cv_single_cell(Y, W, cells, nn_grid, 3,
                                         lambda_unit, lambda_time, tol, max_iter, scores,
                                         "lambda_nn", 1)
        printf("{txt}%s elapsed\n", trop_ms_fmt((trop_ms_now() - t0) / 1000))
        displayflush()
        if (lambda_unit == lu_old & lambda_time == lt_old & lambda_nn == lnn_old) {
            cycles_used = c
            printf("{txt}  converged after %g cycle(s), %s total\n",
                   c, trop_ms_fmt((trop_ms_now() - t0) / 1000))
            displayflush()
            return
        }
    }
    cycles_used = max_cycles
    printf("{txt}  reached max cycles (%g) without full convergence, %s total\n",
           max_cycles, trop_ms_fmt((trop_ms_now() - t0) / 1000))
    displayflush()
}
end

mata:
// Joint grid search on the per-cell criterion. The nn dimension is
// swept warm-started inside trop_loocv_cell_rmse_path
void trop_cv_joint_cell(
    real matrix Y, real matrix W, real matrix cells,
    real rowvector unit_grid, real rowvector time_grid, real rowvector nn_grid,
    real scalar tol, real scalar max_iter,
    real scalar lambda_unit, real scalar lambda_time, real scalar lambda_nn,
    real scalar best_rmse, real scalar n_eval)
{
    real scalar    nu, nt, ng, ii, jj, g, t0
    real scalar    best_u, best_t, best_n
    real rowvector scores, S

    nu = cols(unit_grid); nt = cols(time_grid); ng = cols(nn_grid)
    best_rmse = .; best_u = .; best_t = .; best_n = .
    n_eval = 0
    t0 = trop_ms_now()
    printf("{txt}  joint search: evaluating %g (unit x time) pairs, %g nn values each\n",
           nu * nt, ng)
    displayflush()
    S = trop_ms_init(nu * nt)
    for (ii = 1; ii <= nu; ii++) {
        for (jj = 1; jj <= nt; jj++) {
            scores = trop_loocv_cell_rmse_path(Y, W, cells,
                         unit_grid[ii], time_grid[jj], nn_grid, tol, max_iter)
            n_eval = n_eval + ng
            for (g = 1; g <= ng; g++) {
                if (scores[g] < . & (best_rmse >= . | scores[g] < best_rmse)) {
                    best_rmse = scores[g]
                    best_u = unit_grid[ii]
                    best_t = time_grid[jj]
                    best_n = nn_grid[g]
                }
            }
            trop_ms_tick(S, (ii - 1) * nt + jj, 1)
        }
    }
    if (best_rmse >= .) {
        errprintf("trop_cv_joint_cell(): all CV scores are missing.\n")
        exit(498)
    }
    lambda_unit = best_u; lambda_time = best_t; lambda_nn = best_n
    printf("{txt}  best (%s, %s, %s), %s total\n",
           sprintf("%g", best_u), sprintf("%g", best_t),
           (best_n >= . ? "inf" : sprintf("%g", best_n)),
           trop_ms_fmt((trop_ms_now() - t0) / 1000))
    displayflush()
}
end

mata:
void trop_target_grid()
{
    real matrix    info, TT
    real rowvector taus
    real colvector uvals, tvals
    real scalar    G, U, P, g, ui, tj

    info = st_matrix("__trop_tinfo")        // G x 2, actual unit/time values
    taus = st_matrix("__trop_ttau")         // 1 x G, already rescaled
    G = rows(info)
    uvals = uniqrows(info[., 1])
    tvals = uniqrows(info[., 2])
    U = rows(uvals); P = rows(tvals)
    TT = J(U, P, .)
    for (g = 1; g <= G; g++) {
        ui = selectindex(uvals :== info[g, 1])
        tj = selectindex(tvals :== info[g, 2])
        TT[ui, tj] = taus[g]
    }
    st_matrix("__trop_tgrid", TT)
    st_matrixrowstripe("__trop_tgrid", (J(U, 1, ""), "u" :+ strofreal(uvals)))
    st_matrixcolstripe("__trop_tgrid", (J(P, 1, ""), "t" :+ strofreal(tvals)))
}
end

mata:
string scalar trop_adkey(real scalar u, real scalar t, real scalar n)
{
    return(strofreal(u, "%21.16g") + "|" + strofreal(t, "%21.16g") + "|"
         + strofreal(n, "%21.16g"))
}
end

mata:
string scalar trop_ad_num(real scalar x)
{
    if (x >= .) return("inf")
    return(strtrim(strofreal(x, "%9.3g")))
}
end

mata:
void trop_ad_rangeline(string scalar name,
    real scalar flo, real scalar fhi,
    real scalar alo, real scalar ahi,
    real scalar tick, real scalar best,
    string scalar extra, real scalar ispan)
{
    real scalar W, k, pa, pb, pt, pz
    string rowvector v
    string scalar line

    W = 34
    if (ispan < . & ispan > 0) W = round(34 * (fhi - flo) / ispan)
    if (W < 34) W = 34
    if (W > 68) W = 68
    if (fhi <= flo) return
    v = J(1, W, "-")
    pa = 1 + round((alo - flo) / (fhi - flo) * (W - 1))
    pb = 1 + round((ahi - flo) / (fhi - flo) * (W - 1))
    if (pa < 1)     pa = 1
    if (pa > W - 1) pa = W - 1
    if (pb <= pa)   pb = pa + 1
    if (pb > W) {
        pb = W
        if (pa >= pb) pa = pb - 1
    }
    for (k = pa; k <= pb; k++) v[k] = "="
    v[pa] = "["
    v[pb] = "]"
    if (tick < .) {
        pt = 1 + round((tick - flo) / (fhi - flo) * (W - 1))
        if (pt >= 1 & pt <= W & pt != pa & pt != pb) v[pt] = "|"
    }
    if (best < .) {
        pz = 1 + round((best - flo) / (fhi - flo) * (W - 1))
        if (pz < pa + 1) pz = pa + 1
        if (pz > pb - 1) pz = pb - 1
        if (pz > pa & pz < pb) v[pz] = "*"
    }
    line = ""
    for (k = 1; k <= W; k++) line = line + v[k]
    printf("{txt}    %-12s %s %s %s%s\n",
           name, trop_ad_num(flo), line, trop_ad_num(fhi), extra)
    displayflush()
}
end

mata:
string scalar trop_ad_pad(string scalar s, real scalar w)
{
    string scalar out
    out = s
    while (strlen(out) < w) out = out + " "
    return(out)
}
end

mata:
string scalar trop_ad_axlist(real scalar bd_u, real scalar bd_t, real scalar bd_n)
{
    real scalar k, n
    string rowvector nm
    string scalar out

    nm = J(1, 0, "")
    if (bd_u) nm = nm, "Lambda unit"
    if (bd_t) nm = nm, "Lambda time"
    if (bd_n) nm = nm, "Lambda nn"
    n = cols(nm)
    out = ""
    for (k = 1; k <= n; k++) {
        if (k == 1)      out = nm[k]
        else if (k == n) out = out + " and " + nm[k]
        else             out = out + ", " + nm[k]
    }
    return(out)
}
end

mata:
void trop_ad_heat(transmorphic A,
    real rowvector axis_u, real rowvector axis_t, real rowvector axis_n,
    real scalar bu, real scalar bt,
    real scalar rescale, real scalar gmin, real scalar gmax)
{
    real scalar nu, nt, ng, i, j, g, s, v, lo, hi, idx, bi, bj
    real scalar k, c, pos, prev, col, LW, li, nl, dw, bw, slot, start
    real matrix Sm
    string scalar ramp, line, lab, ch, ticks, left
    string rowvector leg, sym, des

    ramp = "#+-."                     // low RMSE -> "#" ... high RMSE -> "."
    nu = cols(axis_u); nt = cols(axis_t); ng = cols(axis_n)

    Sm = J(nt, nu, .)
    for (j = 1; j <= nt; j++) {
        for (i = 1; i <= nu; i++) {
            v = .
            for (g = 1; g <= ng; g++) {
                s = asarray(A, trop_adkey(axis_u[i], axis_t[j], axis_n[g]))
                if (s < . & (v >= . | s < v)) v = s
            }
            Sm[j, i] = v
        }
    }
    if (rescale) {
        lo = min(Sm); hi = max(Sm)
    }
    else {
        lo = min((gmin, min(Sm))); hi = max((gmax, max(Sm)))
        gmin = lo; gmax = hi
    }
    if (hi <= lo) hi = lo + 1e-12
    bi = 0; bj = 0
    for (i = 1; i <= nu; i++) if (reldif(axis_u[i], bu) < 1e-9) bi = i
    for (j = 1; j <= nt; j++) if (reldif(axis_t[j], bt) < 1e-9) bj = j

    // Legend rows: the ramp maps a cell to round(3*(v-lo)/(hi-lo)), so the
    // symbol boundaries sit at 1/6, 1/2 and 5/6 of the [lo, hi] span. X
    // reports the lowest RMSE on this grid.
    sym = (" ", "X", "#", "+", "-", ".")
    des = ("RMSE",
           "Best RMSE (" + trop_ad_num(min(Sm)) + ")",
           trop_ad_num(lo) + " to " + trop_ad_num(lo + (hi - lo) / 6),
           trop_ad_num(lo + (hi - lo) / 6) + " to " + trop_ad_num(lo + (hi - lo) / 2),
           trop_ad_num(lo + (hi - lo) / 2) + " to " + trop_ad_num(lo + 5 * (hi - lo) / 6),
           trop_ad_num(lo + 5 * (hi - lo) / 6) + " to " + trop_ad_num(hi))
    dw = 0
    for (k = 1; k <= cols(des); k++) dw = max((dw, strlen(des[k])))
    bw = dw + 8
    nl = cols(sym) + 2
    leg = J(1, nl, "")
    line = "+- Legend "
    while (strlen(line) < bw - 1) line = line + "-"
    leg[1] = line + "+"
    for (k = 1; k <= cols(sym); k++)
        leg[k + 1] = "| " + sym[k] + " | " + trop_ad_pad(des[k], dw) + " |"
    line = "+"
    while (strlen(line) < bw - 1) line = line + "-"
    leg[nl] = line + "+"
    li = 1
    LW = 2 * nu + 13                  // width of the heatmap's left column

    // Bottom-align: attachment slots are title + nt rows + border = nt + 2;
    // start attaching so the legend's last line lands on the border.
    slot  = 0
    start = nt + 2 - nl + 1
    if (start < 1) start = 1

    printf("{txt}Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)\n")
    if (rescale)
        printf("{txt}     Each cell shows the lowest RMSE across all lambda_nn's searched in the adapted range.\n")
    else
        printf("{txt}     Each cell shows the lowest RMSE across all lambda_nn's searched.\n")
    printf("{txt}\n")
    left = "Lambda time"
    slot++
    if (slot >= start & li <= nl) {
        printf("{txt}%s      %s\n", trop_ad_pad(left, LW), leg[li])
        li++
    }
    else printf("{txt}%s\n", left)
    for (j = nt; j >= 1; j--) {
        line = ""
        for (i = 1; i <= nu; i++) {
            if (i == bi & j == bj) ch = "X"
            else if (Sm[j, i] >= .) ch = "?"
            else {
                idx = round((Sm[j, i] - lo) / (hi - lo) * 3)
                if (idx < 0) idx = 0
                if (idx > 3) idx = 3
                ch = substr(ramp, idx + 1, 1)
            }
            line = line + " " + ch
        }
        // label the top and bottom rows only
        lab = ((j == nt | j == 1) ? trop_ad_num(axis_t[j]) : "")
        left = sprintf("%9s |%s |", lab, line)
        slot++
        if (slot >= start & li <= nl) {
            printf("{txt}%s      %s\n", trop_ad_pad(left, LW), leg[li])
            li++
        }
        else printf("{txt}%s\n", left)
    }
    line = ""
    for (c = 1; c <= 2 * nu + 1; c++) line = line + "-"
    left = "          +" + line + "+"
    slot++
    if (slot >= start & li <= nl) {
        printf("{txt}%s      %s\n", trop_ad_pad(left, LW), leg[li])
        li++
    }
    else printf("{txt}%s\n", left)

    ticks = ""
    prev = 0
    for (k = 1; k <= 2; k++) {
        if (k == 2 & nu < 2) break
        i = (k == 1 ? 1 : nu)
        lab = trop_ad_num(axis_u[i])
        col = 11 + 2 * i
        pos = col - floor((strlen(lab) + 1) / 2)
        if (pos < prev + 2) pos = prev + 2
        for (c = prev + 1; c < pos; c++) ticks = ticks + " "
        ticks = ticks + lab
        prev = pos + strlen(lab) - 1
    }
    left = ticks + "    Lambda unit"
    if (li <= nl) {
        printf("{txt}%s      %s\n", trop_ad_pad(left, LW), leg[li])
        li++
    }
    else printf("{txt}%s\n", left)
    while (li <= nl) {
        printf("{txt}%s      %s\n", trop_ad_pad("", LW), leg[li])
        li++
    }
    printf("{txt}\n")
    displayflush()
}
end

mata:
void trop_adaptive_eval(
    real matrix Y, real matrix W,
    pointer(real colvector) rowvector sets, real matrix cells, real scalar is_cell,
    real rowvector axis_u, real rowvector axis_t, real rowvector axis_n,
    transmorphic A, real scalar tol, real scalar max_iter,
    real scalar bu, real scalar bt, real scalar bn, real scalar bs)
{
    real scalar    nu, nt, ng, ii, jj, g, allc, s
    real rowvector scores, S
    string scalar  key

    nu = cols(axis_u); nt = cols(axis_t); ng = cols(axis_n)
    S = trop_ms_init(nu * nt)
    for (ii = 1; ii <= nu; ii++) {
        for (jj = 1; jj <= nt; jj++) {
            allc = 1
            for (g = 1; g <= ng; g++) {
                key = trop_adkey(axis_u[ii], axis_t[jj], axis_n[g])
                if (!asarray_contains(A, key)) {
                    allc = 0
                    break
                }
            }
            if (!allc) {
                if (is_cell) {
                    scores = trop_loocv_cell_rmse_path(Y, W, cells,
                                 axis_u[ii], axis_t[jj], axis_n, tol, max_iter)
                }
                else {
                    scores = trop_placebo_rmse_path(Y, sets,
                                 axis_u[ii], axis_t[jj], axis_n, tol, max_iter, 1)
                }
                for (g = 1; g <= ng; g++) {
                    key = trop_adkey(axis_u[ii], axis_t[jj], axis_n[g])
                    asarray(A, key, scores[g])
                }
            }
            trop_ms_tick(S, (ii - 1) * nt + jj, 1)
        }
    }
    bu = .; bt = .; bn = .; bs = .
    for (ii = 1; ii <= nu; ii++) {
        for (jj = 1; jj <= nt; jj++) {
            for (g = 1; g <= ng; g++) {
                s = asarray(A, trop_adkey(axis_u[ii], axis_t[jj], axis_n[g]))
                if (s < . & (bs >= . | s < bs)) {
                    bs = s
                    bu = axis_u[ii]; bt = axis_t[jj]; bn = axis_n[g]
                }
            }
        }
    }
}
end

mata:
// Position-indexed display of a sweep grid, with the winning value bracketed.
void trop_ad_gridline(string scalar name, real rowvector grid, real scalar best)
{
    real scalar   k, hit
    string scalar line, lab

    line = ""
    for (k = 1; k <= cols(grid); k++) {
        lab = trop_ad_num(grid[k])
        if (grid[k] >= .) hit = (best >= .)
        else              hit = (best < . & reldif(grid[k], best) < 1e-9)
        if (hit) line = line + " [" + lab + "]"
        else     line = line + " "  + lab
    }
    printf("{txt}    %-12s%s\n", name, line)
    displayflush()
}
end

mata:
void trop_ad_bracket(real rowvector s, real scalar b, real scalar edge_factor,
                     real scalar lo, real scalar hi, real scalar tick)
{
    real scalar ns, j, k

    tick = .
    ns = cols(s)
    if (ns < 2) {
        lo = s[1]; hi = s[1]
        return
    }
    k = 0
    for (j = 1; j <= ns; j++) {
        if (reldif(s[j], b) < 1e-9) {
            k = j
            break
        }
    }
    if (k == 0) {                        // no match: keep the full span
        lo = s[1]; hi = s[ns]
        return
    }
    if (k == 1) {
        hi = s[2]
        if (s[1] > 0) {
            tick = s[1]
            lo   = max((s[1] - edge_factor * (s[2] - s[1]), 0))
        }
        else lo = s[1]
    }
    else if (k == ns) {
        lo   = s[ns-1]
        tick = s[ns]
        hi   = s[ns] + edge_factor * (s[ns] - s[ns-1])
    }
    else {
        lo = s[k-1]
        hi = s[k+1]
    }
}
end

mata:
real rowvector trop_cv_adaptive(
    real matrix Y, real matrix W,
    pointer(real colvector) rowvector sets, real matrix cells, real scalar is_cell,
    real rowvector g_u, real scalar tune_u,
    real rowvector g_t, real scalar tune_t,
    real rowvector g_n, real scalar tune_n, real scalar nn_inf_only,
    real scalar npts, real scalar max_exp,
    real scalar tol, real scalar max_iter, real scalar vb,
    real scalar phases_used)
{
    real scalar    phase, bu, bt, bn, bs, bd_u, bd_t, bd_n
    real scalar    span, half, nl, nh, zoomed, t0, done, resumed
    real scalar    tick_u, tick_t, tick_n, gmin, gmax
    real scalar    edge_factor, nn_has_inf, fl, fh, k, sc, chg
    real scalar    lo_u, hi_u, lo_t, hi_t, lo_n, hi_n
    real scalar    flo_u, fhi_u, flo_t, fhi_t, flo_n, fhi_n
    real scalar    ispan_u, ispan_t, ispan_n
    real rowvector axis_u, axis_t, axis_n, su, stg, sn
    string colvector keys
    string rowvector parts
    transmorphic   A

    A  = asarray_create()
    edge_factor = 4
    bu = .; bt = .; bn = .; bs = .
    tick_u = .; tick_t = .; tick_n = .
    gmin = .; gmax = .
    phases_used = 0
    t0 = trop_ms_now()
    done = 0

    // Sorted unique finite sweep grids. Fixed lambdas arrive as single-value
    // grids (the command body collapses them before CV starts).
    su  = uniqrows(select(g_u, g_u :< .)')'
    stg = uniqrows(select(g_t, g_t :< .)')'
    nn_has_inf = (sum(g_n :>= .) > 0)
    if (nn_inf_only) sn = J(1, 0, .)
    else             sn = uniqrows(select(g_n, g_n :< .)')'

    ispan_u = (cols(su)  > 1 ? su[cols(su)]   - su[1]  : .)
    ispan_t = (cols(stg) > 1 ? stg[cols(stg)] - stg[1] : .)
    ispan_n = (cols(sn)  > 1 ? sn[cols(sn)]   - sn[1]  : .)

    //--- Phase 1: sweep the supplied grids literally ---//
    axis_u = su
    axis_t = stg
    if (nn_inf_only)     axis_n = (.)
    else if (nn_has_inf) axis_n = sn, (.)
    else                 axis_n = sn

    printf("{txt}Adaptive CV: phase 1 sweeps the %gx%gx%g starting grid of lambdas;\n",
           cols(axis_u), cols(axis_t), cols(axis_n))
    printf("{txt}later phases refine around the winner; points() and expansions() control their cost.\n")
    printf("{txt} \n")
    printf("{txt}  Phase 1: Sweeping the %gx%gx%g (unit, time, nuclear norm) starting grid.\n",
           cols(axis_u), cols(axis_t), cols(axis_n))
    printf("{txt}\n")
    displayflush()
    trop_adaptive_eval(Y, W, sets, cells, is_cell,
                       axis_u, axis_t, axis_n, A, tol, max_iter,
                       bu, bt, bn, bs)
    if (bs >= .) {
        errprintf("trop_cv_adaptive(): all CV scores are missing.\n")
        exit(498)
    }
    phases_used = 1

    bd_u = (tune_u & cols(su) > 1 &
            ((reldif(bu, su[cols(su)]) < 1e-9) |
             (reldif(bu, su[1]) < 1e-9 & su[1] > 0)))
    bd_t = (tune_t & cols(stg) > 1 &
            ((reldif(bt, stg[cols(stg)]) < 1e-9) |
             (reldif(bt, stg[1]) < 1e-9 & stg[1] > 0)))
    bd_n = (tune_n & !nn_inf_only & bn < . & cols(sn) > 1 &
            ((reldif(bn, sn[cols(sn)]) < 1e-9) |
             (reldif(bn, sn[1]) < 1e-9 & sn[1] > 0)))
    printf("{txt}  Phase 1 chose lambdas (%s, %s, %s)%s\n",
           trop_ad_num(bu), trop_ad_num(bt), trop_ad_num(bn),
           ((bd_u | bd_t | bd_n) ? ", on the edge of the swept grid" : ""))
    displayflush()
    if (tune_u) trop_ad_gridline("Lambda unit", axis_u, bu)
    if (tune_t) trop_ad_gridline("Lambda time", axis_t, bt)
    if (tune_n & !nn_inf_only) trop_ad_gridline("Lambda nn", axis_n, bn)
    if (vb) {
        printf("{txt}\n")
        trop_ad_heat(A, axis_u, axis_t, axis_n, bu, bt, 0, gmin, gmax)
    }

    //--- Bracket the winner between its neighbouring grid points ---//
    lo_u = su[1];  hi_u = su[cols(su)]
    lo_t = stg[1]; hi_t = stg[cols(stg)]
    if (nn_inf_only) {
        lo_n = .; hi_n = .
    }
    else {
        lo_n = sn[1]; hi_n = sn[cols(sn)]   // bn = inf keeps the full finite span
    }
    if (tune_u) trop_ad_bracket(su, bu, edge_factor, lo_u, hi_u, tick_u)
    if (tune_t) trop_ad_bracket(stg, bt, edge_factor, lo_t, hi_t, tick_t)
    if (tune_n & !nn_inf_only & bn < .)
        trop_ad_bracket(sn, bn, edge_factor, lo_n, hi_n, tick_n)

    if ((tune_u & cols(su) > 1) | (tune_t & cols(stg) > 1) |
        (tune_n & !nn_inf_only & cols(sn) > 1)) {
        if (!vb) printf("{txt}\n")
        printf("{txt} Adapting search grid:\n")
        if (tune_u & cols(su) > 1) {
            fl = min((su[1], lo_u)); fh = max((su[cols(su)], hi_u))
            trop_ad_rangeline("Lambda unit", fl, fh, lo_u, hi_u, tick_u, .,
                "   [" + trop_ad_num(lo_u) + " ; " + trop_ad_num(hi_u) + "]",
                ispan_u)
        }
        if (tune_t & cols(stg) > 1) {
            fl = min((stg[1], lo_t)); fh = max((stg[cols(stg)], hi_t))
            trop_ad_rangeline("Lambda time", fl, fh, lo_t, hi_t, tick_t, .,
                "   [" + trop_ad_num(lo_t) + " ; " + trop_ad_num(hi_t) + "]",
                ispan_t)
        }
        if (tune_n & !nn_inf_only & cols(sn) > 1) {
            fl = min((sn[1], lo_n)); fh = max((sn[cols(sn)], hi_n))
            trop_ad_rangeline("Lambda nn", fl, fh, lo_n, hi_n, tick_n, .,
                "   [" + trop_ad_num(lo_n) + " ; " + trop_ad_num(hi_n) + "]"
                + (nn_has_inf ? " (and inf)" : ""), ispan_n)
        }
        printf("{txt}\n")
        displayflush()
    }

    while (!done) {

    for (phase = phases_used; phase <= max_exp; phase++) {
        axis_u = (tune_u & hi_u > lo_u ? rangen(lo_u, hi_u, npts)' : (lo_u))
        axis_t = (tune_t & hi_t > lo_t ? rangen(lo_t, hi_t, npts)' : (lo_t))
        if (nn_inf_only) axis_n = (.)
        else {
            axis_n = (tune_n & hi_n > lo_n ? rangen(lo_n, hi_n, npts)' : (lo_n))
            if (tune_n & nn_has_inf) axis_n = axis_n, (.)
        }

        printf("{txt}  Phase %g: Evaluating %gx%gx%g (unit, time, nuclear norm) grid.\n",
               phase + 1, cols(axis_u), cols(axis_t), cols(axis_n))
        printf("{txt}\n")
        displayflush()
        trop_adaptive_eval(Y, W, sets, cells, is_cell,
                           axis_u, axis_t, axis_n, A, tol, max_iter,
                           bu, bt, bn, bs)
        if (bs >= .) {
            errprintf("trop_cv_adaptive(): all CV scores are missing.\n")
            exit(498)
        }

        bd_u = (tune_u & hi_u > lo_u &
                ((reldif(bu, hi_u) < 1e-9) | (reldif(bu, lo_u) < 1e-9 & lo_u > 0)))
        bd_t = (tune_t & hi_t > lo_t &
                ((reldif(bt, hi_t) < 1e-9) | (reldif(bt, lo_t) < 1e-9 & lo_t > 0)))
        bd_n = (tune_n & !nn_inf_only & bn < . & hi_n > lo_n &
                ((reldif(bn, hi_n) < 1e-9) | (reldif(bn, lo_n) < 1e-9 & lo_n > 0)))
        phases_used = phase + 1
        printf("{txt}  Phase %g chose lambdas (%s, %s, %s)%s\n",
               phase + 1,
               trop_ad_num(bu), trop_ad_num(bt), trop_ad_num(bn),
               ((bd_u | bd_t | bd_n) ? ", on the edge of the searched range" : ""))
        displayflush()
        if (tune_u) trop_ad_rangeline("Lambda unit", lo_u, hi_u, lo_u, hi_u, tick_u, bu, "", ispan_u)
        if (tune_t) trop_ad_rangeline("Lambda time", lo_t, hi_t, lo_t, hi_t, tick_t, bt, "", ispan_t)
        if (tune_n & !nn_inf_only)
            trop_ad_rangeline("Lambda nn", lo_n, hi_n, lo_n, hi_n, tick_n, bn,
                              "", ispan_n)
        if (vb) {
            printf("{txt}\n")
            trop_ad_heat(A, axis_u, axis_t, axis_n, bu, bt, 0, gmin, gmax)
        }

        if (!(bd_u | bd_t | bd_n)) break
        if (phase == max_exp) {
            printf("{txt}  Still on the edge of the range after %g expansion(s): the best lambdas may lie outside it.\n",
                   max_exp)
            printf("{txt}  Widen the search with unit_grid(), time_grid() or nn_grid(), or allow more expansions().\n")
            displayflush()
            break
        }
        if (!vb) printf("{txt}\n")
        printf("{txt} Expanding search grid: %s %s on the edge of the range.\n",
               trop_ad_axlist(bd_u, bd_t, bd_n),
               ((bd_u + bd_t + bd_n) > 1 ? "sit" : "sits"))
        printf("{txt}\n")
        displayflush()
        if (bd_u) {
            span = hi_u - lo_u
            if (reldif(bu, lo_u) < 1e-9) {
                tick_u = lo_u
                lo_u = max((lo_u - span, 0))
            }
            else {
                tick_u = hi_u
                hi_u = hi_u + span
            }
        }
        if (bd_t) {
            span = hi_t - lo_t
            if (reldif(bt, lo_t) < 1e-9) {
                tick_t = lo_t
                lo_t = max((lo_t - span, 0))
            }
            else {
                tick_t = hi_t
                hi_t = hi_t + span
            }
        }
        if (bd_n) {
            span = hi_n - lo_n
            if (reldif(bn, lo_n) < 1e-9) {
                tick_n = lo_n
                lo_n = max((lo_n - span, 0))
            }
            else {
                tick_n = hi_n
                hi_n = hi_n + span
            }
        }
    }

    flo_u = lo_u; fhi_u = hi_u
    flo_t = lo_t; fhi_t = hi_t
    flo_n = lo_n; fhi_n = hi_n
    zoomed = 0
    if (tune_u & hi_u > lo_u) {
        half = (hi_u - lo_u) / 4
        nl = max((bu - half, lo_u)); nh = min((bu + half, hi_u))
        if (nl < nh) {
            lo_u = nl; hi_u = nh; zoomed = 1
        }
    }
    if (tune_t & hi_t > lo_t) {
        half = (hi_t - lo_t) / 4
        nl = max((bt - half, lo_t)); nh = min((bt + half, hi_t))
        if (nl < nh) {
            lo_t = nl; hi_t = nh; zoomed = 1
        }
    }
    if (tune_n & !nn_inf_only & bn < . & hi_n > lo_n) {
        half = (hi_n - lo_n) / 4
        nl = max((bn - half, lo_n)); nh = min((bn + half, hi_n))
        if (nl < nh) {
            lo_n = nl; hi_n = nh; zoomed = 1
        }
    }
    if (zoomed) {
        axis_u = (tune_u & hi_u > lo_u ? rangen(lo_u, hi_u, 9)' : (lo_u))
        axis_t = (tune_t & hi_t > lo_t ? rangen(lo_t, hi_t, 9)' : (lo_t))
        if (nn_inf_only) axis_n = (.)
        else {
            axis_n = (tune_n & hi_n > lo_n ? rangen(lo_n, hi_n, 9)' : (lo_n))
            if (tune_n & nn_has_inf) axis_n = axis_n, (.)
        }
        if (!vb) printf("{txt}\n")
        printf("{txt} Adapting search grid:\n")
        if (tune_u) trop_ad_rangeline("Lambda unit", flo_u, fhi_u, lo_u, hi_u, ., .,
                        "   [" + trop_ad_num(lo_u) + " ; " + trop_ad_num(hi_u) + "]",
                        ispan_u)
        if (tune_t) trop_ad_rangeline("Lambda time", flo_t, fhi_t, lo_t, hi_t, ., .,
                        "   [" + trop_ad_num(lo_t) + " ; " + trop_ad_num(hi_t) + "]",
                        ispan_t)
        if (tune_n & !nn_inf_only)
            trop_ad_rangeline("Lambda nn", flo_n, fhi_n, lo_n, hi_n, ., .,
                        "   [" + trop_ad_num(lo_n) + " ; " + trop_ad_num(hi_n) + "]"
                        + (nn_has_inf ? " (and inf)" : ""),
                        ispan_n)
        printf("{txt}\n")
        printf("{txt}  Evaluating %g (unit x time) pairs in range, %g nn values each\n",
               cols(axis_u) * cols(axis_t), cols(axis_n))
        displayflush()
        trop_adaptive_eval(Y, W, sets, cells, is_cell,
                           axis_u, axis_t, axis_n, A, tol, max_iter,
                           bu, bt, bn, bs)
        if (bs >= .) {
            errprintf("trop_cv_adaptive(): all zoom CV scores are missing.\n")
            exit(498)
        }
        printf("{txt}\n")
        printf("{txt}  Best after zoom (%s, %s, %s), %s total\n",
               trop_ad_num(bu), trop_ad_num(bt), trop_ad_num(bn),
               trop_ms_fmt((trop_ms_now() - t0) / 1000))
        displayflush()
        if (tune_u) trop_ad_rangeline("Lambda unit", flo_u, fhi_u, lo_u, hi_u, ., bu, "", ispan_u)
        if (tune_t) trop_ad_rangeline("Lambda time", flo_t, fhi_t, lo_t, hi_t, ., bt, "", ispan_t)
        if (tune_n & !nn_inf_only)
            trop_ad_rangeline("Lambda nn", flo_n, fhi_n, lo_n, hi_n, ., bn,
                              "", ispan_n)
        if (vb) {
            printf("{txt}\n")
            trop_ad_heat(A, axis_u, axis_t, axis_n, bu, bt, 1, gmin, gmax)
        }
    }

    bd_u = (tune_u & fhi_u > flo_u &
            ((reldif(bu, fhi_u) < 1e-9) | (reldif(bu, flo_u) < 1e-9 & flo_u > 0)))
    bd_t = (tune_t & fhi_t > flo_t &
            ((reldif(bt, fhi_t) < 1e-9) | (reldif(bt, flo_t) < 1e-9 & flo_t > 0)))
    bd_n = (tune_n & !nn_inf_only & bn < . & fhi_n > flo_n &
            ((reldif(bn, fhi_n) < 1e-9) | (reldif(bn, flo_n) < 1e-9 & flo_n > 0)))

    if (!(bd_u | bd_t | bd_n)) done = 1
    else if (phases_used > max_exp) {
        printf("{txt}  After zooming, %s %s on the edge of the range and the expansions() budget (%g) is spent.\n",
               trop_ad_axlist(bd_u, bd_t, bd_n),
               ((bd_u + bd_t + bd_n) > 1 ? "sit" : "sits"), max_exp)
        printf("{txt}  The best lambdas may lie outside it; widen the search with unit_grid(), time_grid() or nn_grid().\n")
        displayflush()
        done = 1
    }
    else {
        // Restore the phase-final range (lo_*/hi_* currently hold the zoom
        // window), then expand the offending axis by one span and resume.
        lo_u = flo_u; hi_u = fhi_u
        lo_t = flo_t; hi_t = fhi_t
        lo_n = flo_n; hi_n = fhi_n
        resumed = 0
        if (bd_u) {
            span = hi_u - lo_u
            if (reldif(bu, lo_u) < 1e-9) {
                tick_u = lo_u; lo_u = max((lo_u - span, 0))
            }
            else {
                tick_u = hi_u; hi_u = hi_u + span
            }
            resumed = 1
        }
        if (bd_t) {
            span = hi_t - lo_t
            if (reldif(bt, lo_t) < 1e-9) {
                tick_t = lo_t; lo_t = max((lo_t - span, 0))
            }
            else {
                tick_t = hi_t; hi_t = hi_t + span
            }
            resumed = 1
        }
        if (bd_n) {
            span = hi_n - lo_n
            if (reldif(bn, lo_n) < 1e-9) {
                tick_n = lo_n; lo_n = max((lo_n - span, 0))
            }
            else {
                tick_n = hi_n; hi_n = hi_n + span
            }
            resumed = 1
        }
        if (!resumed) done = 1
        else {
            printf("{txt} Expanding search grid: after zooming, %s %s on the edge of the range.\n",
                   trop_ad_axlist(bd_u, bd_t, bd_n),
                   ((bd_u + bd_t + bd_n) > 1 ? "sit" : "sits"))
            printf("{txt}\n")
            displayflush()
        }
    }

    }   // end outer while

    chg = 0
    keys = asarray_keys(A)
    for (k = 1; k <= rows(keys); k++) {
        sc = asarray(A, keys[k])
        if (sc < . & sc < bs) {
            parts = tokens(subinstr(keys[k], "|", " "))
            bs = sc
            bu = strtoreal(parts[1])
            bt = strtoreal(parts[2])
            bn = strtoreal(parts[3])
            chg = 1
        }
    }
    if (chg) {
        printf("{txt}  Best over all evaluated lambdas (%s, %s, %s)\n",
               trop_ad_num(bu), trop_ad_num(bt), trop_ad_num(bn))
        displayflush()
    }

    return((bu, bt, bn))
}
end
