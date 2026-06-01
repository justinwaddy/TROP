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
trop Y S T D [if] [in] [, group(type) lambda_unit(#) lambda_time(#) lambda_nn(#)
                          cv(method [search] [, suboptions])
                          vce(vcetype [, reps(#) seed(#)]) level(#) verbose ]
```
+ group(): **cell**, or **time**. Determines how treated unit-time cells are grouped into estimands. **cell** (default) treats every unit-time cell as its own target and reports the ATT by aggregating per-cell effects (paper Eq. 1). **time** groups estimands across adjacent treated periods for a given unit, and further groups treatment (blocks) across units which share the same treatment timing.
+ lambda_unit(): tuning parameter for unit weights $\omega_j$. Larger values concentrate weight on control units whose pre-treatment paths most resemble the treated unit; 0 weights all units equally. If omitted, chosen by cross-validation.
+ lambda_time(): tuning parameter for time weights $\theta_s$. Larger values concentrate weight on periods near treatment; 0 weights all periods equally. If omitted, chosen by cross-validation.
+ lambda_nn(): tuning parameter for the nuclear-norm penalty on the low-rank component $\mathbf{L}$. Larger values shrink $\mathbf{L}$ toward low rank; **inf** (or **.**) drops $\mathbf{L}$, reducing TROP to weighted two-way fixed effects. If omitted, chosen by cross-validation.
+ cv(method [search] [, suboptions]): cross-validation scheme used to choose any lambda left unspecified. The first word is the method, the optional second word is the search type, and tuning knobs follow a comma.
  - *method*: **loocv** (default) leave-one-out, holding out each control observation in turn; **resample** random placebo-treated sets drawn from the control panel; **kfold** controls partitioned into folds, each treated once.
  - *search*: **cycle** (default) coordinate descent; **joint** full grid search.
  - *suboptions*:
    - trials(#): placebo draws under **resample** (default 200; ignored otherwise).
    - ntreated(#): placebo-treated units per **resample** draw (default 1); set near your number of actually-treated units so the placebo mimics the real treated-group size.
    - folds(#): number of folds under **kfold** (default 5; ignored otherwise).
    - seed(#): seed under **resample** and **kfold** draws (default 0). Has no effect under loocv, which is deterministic.
    - unit_grid(numlist): candidate $\lambda_\text{unit}$ values. Default `0 0.1 0.2 0.3 0.5 0.8 1.2 1.6 2`.
    - time_grid(numlist): candidate $\lambda_\text{time}$ values. Default `0 0.025 0.05 0.1 0.2 0.35 0.5 0.75 1 2 4`.
    - nn_grid(string): candidate $\lambda_{nn}$ values; may include **.** for the no-$\mathbf{L}$ (TWFE) case. Default `0.005 0.01 0.025 0.05 0.1 0.25 0.5 1 .`.

  E.g. `cv(loocv joint)`, `cv(resample, trials(500) ntreated(10))`, `cv(kfold, folds(10))`, `cv(loocv, nn_grid(0.01 0.1 1 .) unit_grid(0 0.5 1))`.
+ vce(vcetype [, reps(#) seed(#)]): **bootstrap** (default) for stratified block-bootstrap standard errors, or **noinference** to skip inference. reps() sets bootstrap repetitions (default 200); seed() sets the bootstrap seed.
+ level(): confidence level for the reported interval (default 95).
+ verbose: display cross-validation and bootstrap progress.

## Grouping Estimands: Time-block grouping under `group(time)`

By default (`group(cell)`) TROP estimates a separate effect for every treated unit-time cell and averages them. With `group(time)`, treated cells that share the **same contiguous treatment spell** are pooled into a single timing block, one effect is estimated per block, and the blocks are averaged using their treated-cell counts. For example, simultaneous adoption which switches on permanently results in one estimand for $\tau$.

For each unit, find its uninterrupted runs of treated periods. Two runs belong to the same block $g$ if they have identical start and end dates $(a_g, b_g)$. The block is the set of treated cells covered by that timing:

```math
\mathcal{D}_g = \{(i,t) \mid t \in [a_g, b_g],\; W_{it} = 1\}.
```

Each block replaces the cell indicator with a block indicator $I^g_{js} = 1\{(j,s)\in\mathcal{D}_g\}$ and solves

```math
\widehat{\tau}_g(\lambda)
= \arg\min_{\tau_g}\,\min_{\alpha,\beta,L}
\sum_{j,s}
\bigl[(1-W_{js}) + W_{js}I^g_{js}\bigr]\,
\omega^g_j\,\theta^g_s\,
\bigl(Y_{js} - \alpha_j - \beta_s - L_{js} - \tau_g I^g_{js}\bigr)^2
+ \lambda_{nn}\lVert L\rVert_*.
```

The weight $(1-W_{js}) + W_{js}I^g_{js}$ keeps every untreated cell and the treated cells in block $g$, while zero-weighting treated cells from other blocks — so each block is estimated as if it were the only treated block. The grouped ATT is the treated-cell-weighted average of the block effects:

```math
\widehat{\tau}_{\mathrm{time}}
= \sum_{g=1}^{G} \pi_g\,\widehat{\tau}_g,
\qquad
\pi_g = \frac{|\mathcal{D}_g|}{|\mathcal{D}|}.
```

#### Example

```math
W =
\begin{array}{c|cccccc}
 & t{=}1 & t{=}2 & t{=}3 & t{=}4 & t{=}5 & t{=}6 \\
i{=}1 & 0 & 0 & 1 & 1 & 0 & 0 \\
i{=}2 & 0 & 0 & 1 & 1 & 0 & 0 \\
i{=}3 & 0 & 0 & 0 & 1 & 1 & 1 \\
i{=}4 & 0 & 1 & 1 & 0 & 1 & 1 \\
i{=}5 & 0 & 0 & 0 & 0 & 0 & 0
\end{array}
```

- `group(cell)`: 11 separate effects, one per treated cell, averaged equally.
- `group(time)`: contiguous spells are `{3,4}` for units 1 and 2 (identical, so pooled), `{4,5,6}` for unit 3, and `{2,3}` plus `{5,6}` for unit 4 — giving four blocks of sizes 4, 3, 2, 2:

```math
G =
\begin{array}{c|cccccc}
 & t{=}1 & t{=}2 & t{=}3 & t{=}4 & t{=}5 & t{=}6 \\
i{=}1 & 0 & 0 & 1 & 1 & 0 & 0 \\
i{=}2 & 0 & 0 & 1 & 1 & 0 & 0 \\
i{=}3 & 0 & 0 & 0 & 2 & 2 & 2 \\
i{=}4 & 0 & 3 & 3 & 0 & 4 & 4 \\
i{=}5 & 0 & 0 & 0 & 0 & 0 & 0
\end{array}
```

```math
\widehat{\tau}_{\mathrm{time}}
= \frac{4\widehat{\tau}_1 + 3\widehat{\tau}_2 + 2\widehat{\tau}_3 + 2\widehat{\tau}_4}{11}.
```

**In short:** `group(cell)` estimates one effect per treated cell; `group(time)` estimates one effect per shared contiguous treatment spell, then averages by treated-cell count.


## Examples

## References
Abadie, A., Diamond, A., & Hainmueller, J. (2010). [Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program](https://doi.org/10.1198/jasa.2009.ap08746). *Journal of the American Statistical Association*, 105(490), 493–505.

Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021). [Synthetic difference-in-differences](https://doi.org/10.1257/aer.20190159). *American Economic Review*, 111(12), 4088–4118.

Athey, S., Imbens, G., Qu, Z., & Viviano, D. (2025). [Triply robust panel estimators](https://arxiv.org/pdf/2508.21536). arXiv preprint arXiv:2508.21536.
