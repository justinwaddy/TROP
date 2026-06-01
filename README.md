
# TROP: Triply Robust Panel Estimators

This Stata package implements Triply Robust Panel Estimators (TROP) following [Athey et al. (2025)](#references). TROP is a panel-data estimator for treatment effects. It combines a flexible model for potential outcomes, with unit weights which upweight units similar to the treated units, and time weights which upweight the time periods closest to treatment.

This version is currently under development. See our [ReadMe for the future implementations](https://github.com/justinwaddy/TROP/blob/main/README_Future.md) for details on the features which will be implemented in this package.

## Overview

### Setting

TROP operates on an $N \times T$ panel of outcomes $\mathbf{Y}$ and binary treatments $\mathbf{W}$. The estimand is the **average treatment effect on the treated (ATT)**:

$$\tau = \frac{\sum_{i,t} W_{it}\bigl(Y_{it}(1) - Y_{it}(0)\bigr)}{\sum_{i,t} W_{it}}$$

Unlike [SC](https://www.mit.edu/~jhainm/synthpage.html) ([Abadie et al., 2010](#references)) and [SDID](https://github.com/Daniel-Pailanir/sdid) ([Arkhangelsky et al., 2021](#references)), TROP accommodates **general assignment patterns** in which units can move into and out of treatment, and there is no requirement of staggered adoption. Note this package requires simulatenous treatment adoption, but future implementations will implement general assignment patterns.

### The Outcome Model

TROP posits a working model for the untreated potential outcome:

$$Y_{it}(0) = \alpha_i + \beta_t + \mathbf{L}_{it} + \varepsilon_{it}, \qquad \mathbb{E}[\varepsilon_{it} \mid \mathbf{L}] = 0$$

where $\alpha_i$ and $\beta_t$ are unit and time fixed effects, and $\mathbf{L}$ is a low-rank latent component (a factor model). In practice $\mathbf{L}$ is recovered via nuclear-norm-penalized regression.

### The Estimator

TROP jointly estimates fixed effects and the low-rank component by solving a **doubly-weighted nuclear-norm penalized regression** over all control observations. The existing version implements a blocked variant of TROP, under which $\tau$ is estimated across block treatment in $W\_{i,t}$

$$(\hat\tau,\hat\mu,\hat\alpha,\hat\beta,\hat{L}) = \arg\min_{\tau,\mu,\alpha,\beta,L}\ \sum_{j=1}^{N}\sum_{s=1}^{T} \delta_{js}\bigl(Y_{js} - \mu - \alpha_j - \beta_s - L_{js} - \tau\,W_{js}\bigr)^2 + \lambda_{nn}\lVert L\rVert_*,$$

with $\delta_{js} = \omega_j\,\theta_s$ the outer product of unit weights $\omega$ and time weights $\theta$. 

This formulation nests existing estimators as special cases, and these can be requested with arguments in `trop`:
- **DID/TWFE**: $\lambda_{nn} = \infty$, uniform weights $\omega_j = \theta_s = 1$
- **Matrix Completion (MC)**: uniform weights, $\lambda_{nn} < \infty$
- **SC / SDID**: $\lambda_{nn} = \infty$, with SC- or SDID-specific unit and time weights

### Weights

TROP uses two sets of exponentially decaying weights, governed by a tuning triplet $\lambda = (\lambda_\text{time}, \lambda_\text{unit}, \lambda_{nn})$:

**Time weights** discount periods further from the treatment date:

$$\theta_s^{i,t}(\lambda) = \exp\Bigl(-\lambda_\text{time} \cdot |t - s|\Bigr)$$

**Unit weights** downweight control units whose pre-treatment trajectories diverge from the treated unit's:

$$\omega_j^{i,t}(\lambda) = \exp\left(-\lambda_\text{unit} \cdot \left(\frac{\sum_{u \neq t}(1-W_{iu})(1-W_{ju})(Y_{iu}-Y_{ju})^2}{\sum_{u \neq t}(1-W_{iu})(1-W_{ju})}\right)^{1/2}\right)$$


### Tuning via placebo cross-validation
We provide tuning options `cv(resample)` and `cv(kfold)`, which use placebo cross-validation as described in the [TROP Python Tutorial](https://github.com/ostasovskyi/TROP-Estimator/blob/main/notebooks/tutorial.ipynb). The idea is to assign the observed treatment pattern of treated units to a set of control units. The set of control units can either be chosen under `CV(resample)` which draws `ntrials` random control samples, or `CV(k-fold)` which partitions the controls into K folds.

Any lambda left unspecified is chosen by placebo cross-validation, which uses the fact that the estimated effect on a control unit permuted as if it were treated should be near zero. Using only the never-treated units, `trop` assigns the treated time pattern to a set of control units, estimates a placebo $\tau$, and scores a lambda triplet by the root-mean-square of the placebo estimates, 

$$Q(\lambda) = \sqrt{\tfrac{1}{J}\sum_j \hat\tau_j^\text{placebo}(\lambda)^2}$$ 

The placebo sets are drawn once and reused for every lambda on the grid. The sampling method is either `resample` (default), which draws `ntrials` random placebo sets each treating `ntreated` control units, or `kfold` which partitions the controls into `kfold` folds, each treated once. The search is either `cycle` (default), coordinate descent over the three lambdas in turn, or `joint`, an exhaustive grid search. Fixing any subset of lambdas directly via the options leaves only the remaining lambdas to be cross-validated.

## Inputs
+ Y: Outcome variable (numeric)
+ S: Unit variable (numeric or string)
+ T: Time variable (numeric)
+ D: Dummy of treatment, equal to 1 if a unit is treated in a given period, and otherwise 0 (numeric)

The panel must be strongly balanced, with no missing values in Y or D, and D must take only the values 0 and 1 with at least one treated and one control unit.

## Syntax
```s
trop Y S T D [if] [in] [, lambda_unit(#) lambda_time(#) lambda_nn(#)
                          cv(method [search]) ntrials(#) ntreated(#) kfold(#) cv_seed(#)
                          unit_grid(numlist) time_grid(numlist) nn_grid(string)
                          vce(vcetype) reps(#) seed(#) level(#) verbose ]
```
+ lambda_unit(): tuning parameter for unit weights $\omega_j$. Larger values concentrate weight on control units whose pre-treatment paths most resemble the treated unit; 0 weights all units equally. If omitted, chosen by cross-validation.
+ lambda_time(): tuning parameter for time weights $\theta_s$. Larger values concentrate weight on periods near treatment; 0 weights all periods equally. If omitted, chosen by cross-validation.
+ lambda_nn(): tuning parameter for the nuclear-norm penalty on the low-rank component $\mathbf{L}$. Larger values shrink $\mathbf{L}$ toward low rank; **inf** (or **.**) drops $\mathbf{L}$, reducing TROP to weighted two-way fixed effects. If omitted, chosen by cross-validation.
+ cv(method [search]): cross-validation scheme used to choose any lambda left unspecified. The first word is the method, **resample** (default) for random placebo-treated sets drawn from the control panel, or **kfold** for controls partitioned into folds. The optional second word is the search, **cycle** (default) for coordinate descent or **joint** for a full grid search.
+ ntrials(): number of placebo draws under **resample** (default 200; ignored under **kfold**).
+ ntreated(): number of placebo-treated control units per **resample** draw (default 1). Set near your number of actually-treated units so the placebo mimics the real treated-group size.
+ kfold(): number of folds under **kfold** (default 5; ignored under **resample**).
+ cv_seed(): seed for the **resample** and **kfold** draws (default 0).
+ unit_grid(): candidate lambda_unit values. Default `0 0.1 0.2 0.3 0.5 0.8 1.2 1.6 2`.
+ time_grid(): candidate lambda_time values. Default `0 0.025 0.05 0.1 0.2 0.35 0.5 0.75 1 2 4`.
+ nn_grid(): candidate lambda_nn values; may include **.** for the no-L (TWFE) case. Default `0.005 0.01 0.025 0.05 0.1 0.25 0.5 1 .`.
+ vce(): **bootstrap** (default) for stratified block-bootstrap standard errors, or **noinference** to skip inference and report the point estimate only.
+ reps(): bootstrap repetitions (default 200).
+ seed(): bootstrap seed (default 0, which leaves the RNG state unchanged).
+ level(): confidence level for the reported interval (default 95).
+ verbose: display cross-validation and bootstrap progress.

## Examples
The examples use the long-form panels shipped with the package (`trop_panel_penn.csv` and `trop_panel_synthetic.csv`), each with columns `unit time y w`. Load and `xtset` first:
```s
import delimited "trop_panel_penn.csv", varnames(1) clear
xtset unit time
```
 
With all three regularizers off, TROP reduces to textbook DID/TWFE:
```s
trop y unit time w, lambda_unit(0) lambda_time(0) lambda_nn(inf) vce(noinference)
```
which returns
```
----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.22172
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          29
-------------+--------------------------------------------------
 lambda_unit |      0.0000
 lambda_time |      0.0000
   lambda_nn |         inf
----------------------------------------------------------------
```
 
Turning on both weight dimensions gives weighted two-way fixed effects:
```s
trop y unit time w, lambda_unit(0.3) lambda_time(0.325) lambda_nn(0.016) vce(noinference)
```
which returns
```
----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.04573
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          29
-------------+--------------------------------------------------
 lambda_unit |      0.3000
 lambda_time |      0.3250
   lambda_nn |        .016
----------------------------------------------------------------

```
 
No weights with a finite nuclear penalty i.e. matrix completion:
```s
trop y unit time w, lambda_unit(0) lambda_time(0) lambda_nn(0.6) vce(noinference)
```
which returns
```
----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.06312
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          29
-------------+--------------------------------------------------
 lambda_unit |      0.0000
 lambda_time |      0.0000
   lambda_nn |          .6
----------------------------------------------------------------
```

## References
Abadie, A., Diamond, A., & Hainmueller, J. (2010). [Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program](https://doi.org/10.1198/jasa.2009.ap08746). *Journal of the American Statistical Association*, 105(490), 493-505.

Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021). [Synthetic difference-in-differences](https://doi.org/10.1257/aer.20190159). *American Economic Review*, 111(12), 4088-4118.

Athey, S., Imbens, G., Qu, Z., & Viviano, D. (2025). [Triply robust panel estimators](https://arxiv.org/pdf/2508.21536). arXiv preprint arXiv:2508.21536.
