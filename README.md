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

where $I_{j,s}^{i,t}$ is equal to 1 when $(j,s) = (i,t)$, and 0 otherwise. For a given target cell $(i,t)$, the estimator recovers its treatment effect as if it were the only treated unit. In other words, the front term $[(1-W_{j,s}) + W_{j,s} I^{i,t}_{j,s}]$ zero-weights every other treated cell, and $\tau\_{i,t}\cdot I^{i,t}\_{j,s}$ recovers the treatment effect for only the given target unit. This generalises the estimator to treated and _untreated (placebo) units_. 

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

over a grid, cycling through the three parameters in sequence. We also provide alternative tuning options `cv(resample)` and `cv(kfold)`, which use placebo cross-validation as described in the [TROP Python Tutorial](https://github.com/ostasovskyi/TROP-Estimator/blob/main/notebooks/tutorial.ipynb). The idea is to assign the observed treatment pattern of treated units to a set of control units. The set of control units can either be chosen under `CV(resample)' which draws `ntrials' random control samples, or CV(k-fold) which partitions the controls into K folds.

## Syntax
```s
trop Y S T D [if] [in], vce(method) seed(#)  
```
where $Y$ is an outcome of interest, $S$ and $T$ are variables indicating unit and time indicators respectively, and $D$ is a binary variable indicating treatment assignment which switches to 1 for units and time periods when treatment is assigned.  

+ vce(): CHECK IF REQURED
+ seed(): seed define for pseudo-random numbers.
...

## Examples


## References
Abadie, A., Diamond, A., & Hainmueller, J. (2010). [Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program](https://doi.org/10.1198/jasa.2009.ap08746). *Journal of the American Statistical Association*, 105(490), 493–505.

Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021). [Synthetic difference-in-differences](https://doi.org/10.1257/aer.20190159). *American Economic Review*, 111(12), 4088–4118.

Athey, S., Imbens, G., Qu, Z., & Viviano, D. (2025). [Triply robust panel estimators](https://arxiv.org/pdf/2508.21536). arXiv preprint arXiv:2508.21536.
