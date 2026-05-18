# TROP: Triply Robust Panel Estimators

This Stata package implements Triply Robust Panel Estimators (TROP) following Athey et al. (2025). TROP is a panel-data estimator for treatment effects. It combines a flexible model for potential outcomes, with unit weights which upweight units similar to the treated units, and time weights which upweight the time periods closest to treatment.

## Overview

### Setting

TROP operates on an $N \times T$ panel of outcomes $\mathbf{Y}$ and binary treatments $\mathbf{W}$. The estimand is the **average treatment effect on the treated (ATT)**:

$$\tau = \frac{\sum_{i,t} W_{it}\bigl(Y_{it}(1) - Y_{it}(0)\bigr)}{\sum_{i,t} W_{it}}$$

Unlike SC and SDID, TROP accommodates **general assignment patterns** — units can move into and out of treatment, and there is no requirement of staggered adoption.

### The Outcome Model

TROP posits a working model for the untreated potential outcome:

$$Y_{it}(0) = \alpha_i + \beta_t + \mathbf{L}_{it} + \varepsilon_{it}, \qquad \mathbb{E}[\varepsilon_{it} \mid \mathbf{L}] = 0$$

where $\alpha_i$ and $\beta_t$ are unit and time fixed effects, and $\mathbf{L}$ is a low-rank latent component (a factor model). In practice $\mathbf{L}$ is recovered via nuclear-norm-penalized regression.

### The Estimator

For a treated unit–period pair $(i^*, t^*)$, TROP jointly estimates fixed effects and the low-rank component by solving a **doubly-weighted nuclear-norm penalized regression** over all control observations:

$$(\hat\alpha, \hat\beta, \hat{\mathbf{L}}) = \operatorname*{arg\,min}_{\alpha,\beta,\mathbf{L}} \sum_{j,s} \theta_s^{i,t}\,\omega_j^{i,t}(1 - W_{js})\bigl(Y_{js} - \alpha_j - \beta_s - L_{js}\bigr)^2 + \lambda_{nn}\|\mathbf{L}\|_*$$

The treatment effect estimate is then:

$$\hat\tau_{it} = Y_{it} - \hat\alpha_i - \hat\beta_t - \hat{\mathbf{L}}_{it}$$

This formulation nests existing estimators as special cases, and these can be requested with arguments in `trop`:
- **DID/TWFE**: $\lambda_{nn} = \infty$, uniform weights $\omega_j = \theta_s = 1$
- **Matrix Completion (MC)**: uniform weights, $\lambda_{nn} < \infty$
- **SC / SDID**: $\lambda_{nn} = \infty$, with SC- or SDID-specific unit and time weights

### Weights

TROP uses two sets of exponentially decaying weights, governed by a tuning triplet $\lambda = (\lambda_\text{time}, \lambda_\text{unit}, \lambda_{nn})$:

**Time weights** discount periods further from the treatment date:

$$\theta_s^{i,t}(\lambda) = \exp\!\Bigl(-\lambda_\text{time} \cdot |t - s|\Bigr)$$

**Unit weights** downweight control units whose pre-treatment trajectories diverge from the treated unit's:

$$\omega_j^{i,t}(\lambda) = \exp\!\left(-\lambda_\text{unit} \cdot \left(\frac{\sum_{u \neq t}(1-W_{iu})(1-W_{ju})(Y_{iu}-Y_{ju})^2}{\sum_{u \neq t}(1-W_{iu})(1-W_{ju})}\right)^{1/2}\right)$$

### Tuning via Leave-One-Out Cross-Validation

The triplet $(\lambda_\text{time}, \lambda_\text{unit}, \lambda_{nn})$ is chosen by a **leave-one-out cross-validation (LOOCV)** criterion that exploits the fact that the estimated treatment effect on any control unit---permuted as if it were the true treated unit---should be close to zero. Formally, TROP minimizes:

$$Q(\lambda) = \sum_{i,t} (1 - W_{it})\,\bigl(\hat\tau_{it}^\text{loocv}(\lambda)\bigr)^2$$

over a grid, cycling through the three parameters in sequence. 

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
Athey, S., Imbens, G., Qu, Z., & Viviano, D. (2025). [Triply robust panel estimators](https://arxiv.org/pdf/2508.21536). arXiv preprint arXiv:2508.21536.
