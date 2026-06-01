# TROP: Triply Robust Panel Estimators

This Stata package implements Triply Robust Panel Estimators (TROP) following [Athey et al. (2025)](#references). TROP is a panel-data estimator for treatment effects. It combines a flexible model for potential outcomes, with unit weights which upweight units similar to the treated units, and time weights which upweight the time periods closest to treatment.

## Overview

### Setting

TROP operates on an $N \times T$ panel of outcomes $\mathbf{Y}$ and binary treatments $\mathbf{W}$. The estimand is the **average treatment effect on the treated (ATT)**:

$$\tau = \frac{\sum_{i,t} W_{it}\bigl(Y_{it}(1) - Y_{it}(0)\bigr)}{\sum_{i,t} W_{it}}$$

Unlike [SC](https://www.mit.edu/~jhainm/synthpage.html) ([Abadie et al., 2010](#references)) and [SDID](https://github.com/Daniel-Pailanir/sdid) ([Arkhangelsky et al., 2021](#references)), TROP accommodates **general assignment patterns** in which units can move into and out of treatment, and there is no requirement of staggered adoption.

### The Outcome Model

TROP posits a working model for the untreated potential outcome:

$$Y_{it}(0) = \alpha_i + \beta_t + \mathbf{L}_{it} + \varepsilon_{it}, \qquad \mathbb{E}[\varepsilon_{it} \mid \mathbf{L}] = 0$$

where $\alpha_i$ and $\beta_t$ are unit and time fixed effects, and $\mathbf{L}$ is a low-rank latent component (a factor model). In practice $\mathbf{L}$ is recovered via nuclear-norm-penalized regression.

### The Estimator

For a treated unit–period pair $(i^{\*}, t^{\*})$, TROP jointly estimates fixed effects and the low-rank component by solving a **doubly-weighted nuclear-norm penalized regression** over all control observations. For each unit $(i,t)$, TROP constructs the estimator

$$\hat\tau_{it}(\lambda) = \underset{\tau_{i,t}}{\arg\min}\ \min_{\alpha,\beta,L}\ \sum_{j=1}^{N}\sum_{s=1}^{T} \bigl[(1-W_{j,s}) + W_{j,s}\,I^{i,t}_{j,s}\bigr]\\omega^{i,t}_j(\lambda)\\theta^{i,t}_s(\lambda)\\bigl(Y_{js}-\alpha_j-\beta_s-L_{js}-\tau_{i,t}\cdot I^{i,t}_{j,s}\bigr)^2 + \lambda_{nn}\lVert L\rVert$$

where $I_{j,s}^{i,t}$ is equal to 1 when $(j,s) = (i,t)$, and 0 otherwise. For a given target cell $(i,t)$, the estimator recovers its treatment effect as if it were the only treated unit. In other words, the front term $[(1-W_{j,s}) + W_{j,s} I^{i,t}_{j,s}]$ zero-weights every other treated cell, and $\tau\_{i,t}\cdot I^{i,t}\_{j,s}$ recovers the treatment effect for only the given target unit. The estimator is specified such that it generalises TROP for treated units and _untreated (placebo) units_. 

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

### Tuning via Leave-One-Out Cross-Validation

By default, the triplet $(\lambda_\text{time}, \lambda_\text{unit}, \lambda_{nn})$ is chosen by a **leave-one-out cross-validation (LOOCV)** criterion that exploits the fact that the estimated treatment effect on any control unit---permuted as if it were the true treated unit---should be close to zero. Formally, TROP minimizes:

$$Q(\lambda) = \sum_{i,t} (1 - W_{it})\,\bigl(\hat\tau_{it}^\text{loocv}(\lambda)\bigr)^2$$

over a grid, cycling through the three parameters in sequence. Users can also specify lambda values directly with `lambda_unit()`, `lambda_time()`, and `lambda_nn()`. Any $λ$ left unspecified is chosen by $LOOCV$ by default.

We also provide alternative tuning options `cv(resample)` and `cv(kfold)`, which use placebo cross-validation as described in the [TROP Python Tutorial](https://github.com/ostasovskyi/TROP-Estimator/blob/main/notebooks/tutorial.ipynb). The idea is to assign the observed treatment pattern of treated units to a set of control units. The set of control units can either be chosen under `CV(resample)` which draws `ntrials` random control samples, or `CV(k-fold)` which partitions the controls into K folds.

## Inputs
+ Y: Outcome variable (numeric)
+ S: Unit variable (numeric or string)
+ T: Time variable (numeric)
+ D: Dummy of treatment, equal to 1 if a unit is treated in a given period, and otherwise 0 (numeric)

The panel must be balanced, with no missing values in Y or D. Treatment may follow any 0/1 pattern: units can move into and out of treatment, with no requirement of staggered or absorbing adoption.

## Syntax
```s
trop Y S T D [if] [in], group(type) lambda_unit(#) lambda_time(#) lambda_nn(#)
                        cv(method [search]) unit_grid(numlist) time_grid(numlist)
                        nn_grid(string) ntrials(#) kfold(#) ntreated(#) cv_seed(#)
                        vce(method) reps(#) seed(#) level(#) verbose
```
+ group(): **cell**, or **time**. Sets how treated unit-time cells are grouped into estimands. **cell** (default) treats every unit-time cell as its own target and reports the cell-count-weighted average of the per-cell effects (paper Eq. 1). **time** groups estimands across cells which have adjacent treated time periods, and further groups treatment (blocks) across units which share the same treatment timing. 
+ vce(): **bootstrap** (default) for stratified block-bootstrap standard errors. If you want to omit this procedure use **noinference**.
+ lambda_unit(): tuning parameter for unit weights $\omega_j$. Larger values concentrate weight on control units whose pre-treatment paths most resemble the treated unit; 0 weights all units equally. If omitted, chosen by cross-validation.
+ lambda_time(): tuning parameter for time weights $\theta_s$. Larger values concentrate weight on periods near treatment; 0 weights all periods equally. If omitted, chosen by cross-validation.
+ lambda_nn(): tuning parameter for nuclear-norm penalty on the low-rank component $\mathbf{L}$. Larger values shrink $\mathbf{L}$ toward low rank; **inf** (or **.**) drops $\mathbf{L}$, reducing TROP to weighted two-way fixed effects. If omitted, chosen by cross-validation.
+ cv(): cross-validation scheme. The first word is the method, the second is the search type. For the method, options incude **loocv** (default) for the cell-level leave-one-out criterion, or **resample**/**kfold** for placebo CV on the control panel. The optional second word is the search type: **cycle** (default) for coordinate descent, or **joint** for a full grid search. E.g. cv(loocv joint), cv(resample).
+ unit_grid(): grid of candidate lambda_unit values searched during cross-validation. Default is 0 0.1 0.2 0.3 0.5 0.8 1.2 1.6 2.
+ time_grid(): grid of candidate lambda_time values searched during cross-validation. Default is 0 0.025 0.05 0.1 0.2 0.35 0.5 0.75 1 2 4.
+ nn_grid(): grid of candidate lambda_nn values searched during cross-validation; may include . for the no-L (TWFE) case. Default is 0.005 0.01 0.025 0.05 0.1 0.25 0.5 1 .
+ ntrials(): number of cross-validation placebo draws (default 200): held-out $W=0$ cells under cv(loocv), or placebo samples under cv(resample).
+ kfold(): number of folds for cv(kfold) (default 5).
+ ntreated(): number of control units assigned placebo treatment in each cv(resample) draw (default 1). Ignored unless cv(resample) is used; set it near your number of actually-treated units so the placebo mimics the real treated-group size.
+ cv_seed(): seed for the cross-validation placebo draws (default 0), for reproducibility.
+ reps(): repetitions for the bootstrap standard errors (default 200).
+ seed(): seed for the bootstrap draws.
+ level(): confidence level for the reported interval (default 95).
+ verbose: display cross-validation and bootstrap progress.

## Examples


## References
Abadie, A., Diamond, A., & Hainmueller, J. (2010). [Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program](https://doi.org/10.1198/jasa.2009.ap08746). *Journal of the American Statistical Association*, 105(490), 493–505.

Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021). [Synthetic difference-in-differences](https://doi.org/10.1257/aer.20190159). *American Economic Review*, 111(12), 4088–4118.

Athey, S., Imbens, G., Qu, Z., & Viviano, D. (2025). [Triply robust panel estimators](https://arxiv.org/pdf/2508.21536). arXiv preprint arXiv:2508.21536.
