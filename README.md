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

Future releases will include arguments to request:
- **SC / SDID**: $\lambda_{nn} = \infty$, with SC- or SDID-specific unit and time weights

### Weights

TROP uses two sets of exponentially decaying weights, governed by a tuning triplet $\lambda = (\lambda_\text{time}, \lambda_\text{unit}, \lambda_{nn})$:

**Time weights** discount periods further from the treatment date:

$$\theta_s^{i,t}(\lambda) = \exp\Bigl(-\lambda_\text{time} \cdot |t - s|\Bigr)$$

**Unit weights** downweight control units whose pre-treatment trajectories diverge from the treated unit's:

$$\omega_j^{i,t}(\lambda) = \exp\left(-\lambda_\text{unit} \cdot \left(\frac{\sum_{u \neq t}(1-W_{iu})(1-W_{ju})(Y_{iu}-Y_{ju})^2}{\sum_{u \neq t}(1-W_{iu})(1-W_{ju})}\right)^{1/2}\right)$$

### Tuning via cross-validation

Users can specify lambda values directly with `lambda_unit()`, `lambda_time()`, and `lambda_nn()`. Any $\lambda$ left unspecified is chosen by cross-validation. CV repeatedly estimates a placebo treatment effect for control units under different combinations of lambdas, to find the combination of lambdas which minimises placebo average treatment effects.

+ Under `group(time)`, the placebo-resampling methods `cv(resample)` (default) and `cv(kfold)` use placebo cross-validation. For `cv(resample)`, we use never-treated control units, assign the treatment pattern of treated units to a random subset, and estimate the average treatment effect per panel. `ntrials` specifies the number of trials. We also provide a `cv(kfold)` option which partitions the controls into `folds=5` (default) folds, assigns treatment to the entire fold (recycling treated unit patterns), and uses held-out units from other folds as the donor pool. The score for resample and k-fold is $$Q(\lambda) = \frac{1}{B}\sum_{b=1}^{B} \bigl(\hat\tau^{\text{placebo}}_{b}(\lambda)\bigr)^2$$ where B is the number of trials or folds for each CV method respectively.


+ Under `group(cell)`, the LOOCV criterion for `cv(loocv)` exploits the fact that the estimated treatment effect on any control unit---permuted as if it were the true treated unit---should be close to zero. Formally, it estimates a treatment effect for each control unit-cell, and minimizes
$$Q(\lambda) = \sum_{i,t} (1 - W_{it})\,\bigl(\hat\tau_{it}^\text{loocv}(\lambda)\bigr)^2$$
over the lambda. 

> **Note** The outcome is standardized ($(Y-\text{mean})/\text{SD}$) before fitting, and $\tau$, the standard error, and the confidence interval are mapped back to the raw outcome scale. The $\lambda$ grids and the returned `lambda_unit`/`lambda_time`/`lambda_nn` are therefore on a standardized-outcome scale. 


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
                          vce(vcetype [, reps(#) seed(#)])
                          level(#) verbose detail ]
```
+ group(): **time** (default) or **cell**. Determines how treated unit-time cells are grouped into estimands. **time** (default) groups estimands across adjacent treated periods for a given unit into blocks, and further groups these blocks across units which share the same treatment timing. For example, simultaneous adoption which switches on permanently results in one estimand for $\tau$, and staggered adoption results in as many treatment effects as there are adoption cohorts as usual. **cell** treats every unit-time cell as its own target and reports the ATT by aggregating per-cell effects (paper Eq. 1).
+ lambda_unit(): tuning parameter for unit weights $\omega_j$. Larger values concentrate weight on control units whose pre-treatment paths most resemble the treated unit; 0 weights all units equally. If omitted, chosen by cross-validation.
+ lambda_time(): tuning parameter for time weights $\theta_s$. Larger values concentrate weight on periods near treatment; 0 weights all periods equally. If omitted, chosen by cross-validation.
+ lambda_nn(): tuning parameter for the nuclear-norm penalty on the low-rank component $\mathbf{L}$. Larger values shrink $\mathbf{L}$ toward low rank; **inf** (or **.**) drops $\mathbf{L}$, reducing TROP to weighted two-way fixed effects. If omitted, chosen by cross-validation.
+ cv(method [search] [, suboptions]): cross-validation scheme used to choose any lambda left unspecified. The first word is the method, the optional second word is the search type, and tuning knobs follow a comma.
  - *method*: **loocv** leave-one-out (default under `group(cell)`); **resample** random placebo-treated sets drawn from the control panel (default under `group(time)`); **kfold** controls partitioned into folds, each treated once. **resample** and **kfold** require `group(time)`, **loocv** requires `group(cell)`.
  - *search*: **cycle** (default) coordinate descent, optimising one lambda at a time over its grid with the other two held at their current values, repeating until the triplet stops changing (at most 50 cycles); **joint** an exhaustive grid search over every combination; **adaptive** evaluates an evenly spaced grid across a range, widens the range whenever the best lambda lands on its edge, then zooms in around the best point and searches again. Uses the smallest and largest values of any grid you supply as its starting range ([0;1] default grid).
  - *suboptions*:
    - trials(#): placebo draws under **resample** (default 200; ignored otherwise).
    - folds(#): number of folds under **kfold** (default 5; ignored otherwise).
    - cells(#): number of randomly sampled control cells scored under **loocv** (default all control cells; ignored otherwise). Useful for large panels where full LOOCV is slow.
    - seed(#): seed for the resample and kfold draws (default 0). Under LOOCV it has an effect only with the cells(#) option. LOOCV with all control cells is deterministic.
    - points(#): grid points per lambda under **adaptive** (default 10).
    - expansions(#): how many times **adaptive** may widen a range whose best lambda sits on edge of grid (default 6).
    - unit_grid(numlist): candidate $\lambda_\text{unit}$ values. Default `0 0.1 0.2 0.3 0.5 0.8 1.2 1.6 2`.
    - time_grid(numlist): candidate $\lambda_\text{time}$ values. Default `0 0.025 0.05 0.1 0.2 0.35 0.5 0.75 1 2 4`.
    - nn_grid(string): candidate $\lambda_{nn}$ values; may include **.** for the no-$\mathbf{L}$ (TWFE) case. Default `0.005 0.01 0.025 0.05 0.1 0.25 0.5 1 .`.

  E.g. `cv(loocv joint)`, `cv(resample, trials(500))`, `cv(kfold, folds(10))`, `cv(resample adaptive, points(12) time_grid(0 4))`, `cv(loocv, nn_grid(0.01 0.1 1 .) unit_grid(0 0.5 1))`.
+ vce(vcetype [, reps(#) seed(#)]): **noinference** (default) to skip inference; **bootstrap** for stratified block-bootstrap standard errors and a percentile confidence interval; **jackknife** for delete-one-unit standard errors. reps() sets bootstrap repetitions (default 200) and seed() sets the bootstrap seed.
+ level(): confidence level for the reported interval (default 95).
+ verbose: under `cv(..., adaptive)`, print an RMSE heat map of the searched grid after each phase, and print the selected lambdas at the end of cross-validation.
+ detail: display which units are grouped into each cohort under `group(time)`.

## Grouping Estimands: Per-cell vs Time

There are two options to group estimands:

+ Under `group(cell)` TROP estimates a separate effect for every treated unit-time cell and averages them. This is faithful to the paper and provides heterogeneous treatment effects, but is computationally expensive relative to grouping estimands.

+ Under `group(time)` (default), treated cells that share the **same uninterrupted treatment period** are pooled into a block across time, blocks which share the same uninterrupted treatment period are further pooled across units, and the treatment effect is averaged. For example, simultaneous adoption which switches on permanently results in one estimand for $\tau$, whereas staggered adoption results in as many estimands as there are cohorts as usual. Time distances are determined using the mid-point of the block. Unit distances are determined relative to the average outcome at each time period of the treated block. For further details, see the [interactive TROP example](https://justinwaddy.co.uk/tropinteractive.html). 

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
- `group(time)`: treatment periods are `{3,4}` for units 1 and 2 (identical, so pooled), `{4,5,6}` for unit 3, and `{2,3}` plus `{5,6}` for unit 4 — giving four blocks of sizes 4, 3, 2, 2:

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

## Examples

In these examples, we use Penn Word Tables (`trop_panel_penn.csv`) to conduct a null effect study. First, load in the data:
```s
import delimited "trop_panel_penn.csv", varnames(1) clear
xtset unit time
```
Next, we will generate four different treatment patterns. The patterns are: one treated unit with a single treated period, simultaneous adoption, staggered adoption, and general pattern assignment.

```s
gen w_single = (unit >= 111) & (time == 40) //Only one treated time period and one treated unit

gen w_block = (unit >= 97) & (time >= 31) //Multiple treated periods and units

gen w_stag = 0
replace w_stag = 1 if inrange(unit,  97, 101) & time >= 21 //Three different treatment cohorts
replace w_stag = 1 if inrange(unit, 102, 106) & time >= 31
replace w_stag = 1 if inrange(unit, 107, 111) & time >= 41

set seed 42
gen w_gen = (unit >= 102) & (runiform() < 0.2) //We choose the same treated units as before, and given treated-unit cell has a 20% chance of being assigned treatment status.
```

### Single treated period and unit
 
First, we estimate the effects for one treated unit treated in the final period. Grouping by cell and time will be equivalent here under the same lambdas. Firstly, we'll test out the command using choice lambdas (which computes almost instantly):

```s
trop y unit time w_single, lambda_unit(0.3) lambda_time(0.325) lambda_nn(0.1)
```

```
----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.02629
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |           1
-------------+--------------------------------------------------
 lambda_unit |      0.3000
 lambda_time |      0.3250
   lambda_nn |          .1
----------------------------------------------------------------
```

No unit is actually treated in this null effect study, so the true effect is zero and the estimate here is small. Using set lambdas are good for a first run, but they were chosen arbitrarily so let's use cross-validation to choose them. We will begin by using resampling cross-validation:
```s
trop y unit time w_single, cv(resample, seed(1))
```
 
```
Cross-validating lambdas using resample with cycle search, and 200 trials (seed 1).
To reduce resample computational time, reduce no of trials or set lambdas.
  cycle 1 of up to 50:  lambda_unit -> .2   lambda_time -> 2   lambda_nn -> .05   4:26 elapsed
  cycle 2 of up to 50:  lambda_unit -> .8   lambda_time -> 2   lambda_nn -> .05   12:11 elapsed
  cycle 3 of up to 50:  lambda_unit -> .8   lambda_time -> 2   lambda_nn -> .05   21:11 elapsed
  converged after 3 cycle(s), 21:11 total

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.02356
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |           1
-------------+--------------------------------------------------
 lambda_unit |      0.8000
 lambda_time |      2.0000
   lambda_nn |         .05
             |  (selected by resample CV)
----------------------------------------------------------------
```
 
This takes about twenty minutes and moves the estimate slightly closer to zero. The cycle search converged after three passes at `lambda_unit = 0.8` and `lambda_time = 2`, both interior to the default grids. With a single treated cell, `group(cell)` and `group(time)` target the same estimand, so LOOCV gives a directly comparable estimate:
```s
trop y unit time w_single, group(cell) cv(loocv, cells(200) seed(1))
```
 
```
Cross-validating lambdas using loocv with cycle search, and 200 samples of 5327 total control cells.
To reduce loocv computational time, reduce number of cells or set lambdas.
  marginal:  lambda_time -> 4   lambda_unit -> 2   lambda_nn -> .25   1:10 elapsed
  cycle 1 of up to 50:  lambda_unit -> 2   lambda_time -> 4   lambda_nn -> .025   4:48 elapsed
  cycle 2 of up to 50:  lambda_unit -> 2   lambda_time -> 4   lambda_nn -> .025   14:46 elapsed
  converged after 2 cycle(s), 14:46 total

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.01512
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |           1
-------------+--------------------------------------------------
 lambda_unit |      2.0000
 lambda_time |      4.0000
   lambda_nn |        .025
             |  (selected by loocv CV)
----------------------------------------------------------------
```

LOOCV lands closer to zero (0.015 against 0.024). LOOCV finished in 14:46 but scored only 200 of the 5,327 control cells, set by `cells(200)`. Full LOOCV on this panel is considerably slower than the resample run. `lambda_unit = 2` and `lambda_time = 4` are the largest values in the default grids, and when cross-validation selects a boundary the grid is probably too narrow, so the reported optimum may be a result of lying on the edge of the grid (rather than at a minimum). 

The `adaptive` search instead sweeps a coarse grid, widens it whenever the winner lands on an edge, and zooms in on the best point. With `verbose` it prints an RMSE heat map after each phase which you can see by expanding the output below:

```s
trop y unit time w_single, cv(resample adaptive, trials(200)) verbose
```

<details>
<summary><b>Full adaptive cross-validation log</b> (4 phases, 6 hours 51 minutes total)</summary>

```
Cross-validating lambdas using resample with adaptive search, and 200 trials (seed 0).
To reduce resample computational time, reduce no of trials or set lambdas.
Adaptive CV: phase 1 sweeps the 9x11x9 starting grid of lambdas;
later phases refine around the winner; points() and expansions() control their cost.
 
  Phase 1: Sweeping the 9x11x9 (unit, time, nuclear norm) starting grid.

  First replication 0:23, roughly 37:57 expected in total
    5%  (5/99)   7:06 elapsed   about 2:13:28 remaining
   10%  (10/99)   18:05 elapsed   about 2:40:56 remaining
   15%  (15/99)   24:44 elapsed   about 2:18:30 remaining
   20%  (20/99)   38:03 elapsed   about 2:30:17 remaining
   25%  (25/99)   44:36 elapsed   about 2:12:00 remaining
   30%  (30/99)   57:44 elapsed   about 2:12:47 remaining
   35%  (35/99)   1:05:04 elapsed   about 1:58:58 remaining
   40%  (40/99)   1:19:08 elapsed   about 1:56:43 remaining
   45%  (45/99)   1:26:15 elapsed   about 1:43:30 remaining
   50%  (50/99)   1:36:49 elapsed   about 1:34:52 remaining
   55%  (55/99)   1:43:22 elapsed   about 1:22:41 remaining
   60%  (60/99)   1:52:36 elapsed   about 1:13:11 remaining
   65%  (65/99)   1:58:43 elapsed   about 1:02:05 remaining
   70%  (70/99)   2:04:33 elapsed   about 51:35 remaining
   75%  (75/99)   2:12:31 elapsed   about 42:24 remaining
   80%  (80/99)   2:18:13 elapsed   about 32:49 remaining
   85%  (85/99)   2:33:07 elapsed   about 25:13 remaining
   90%  (90/99)   2:40:46 elapsed   about 16:04 remaining
   95%  (95/99)   2:54:26 elapsed   about 7:20 remaining
  100%  (99/99)   2:56:38 elapsed
  Phase 1 chose lambdas (1.2, 4, .005), on the edge of the swept grid
    Lambda unit  0 .1 .2 .3 .5 .8 [1.2] 1.6 2
    Lambda time  0 .025 .05 .1 .2 .35 .5 .75 1 2 [4]
    Lambda nn    [.005] .01 .025 .05 .1 .25 .5 1 inf

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
        4 | + # # # # # X # # |
          | + # # # # # # # # |
          | + + + + + + + + # |
          | + + + + + + + # # |
          | + + + + + + + + - |      +- Legend --------------+
          | + + + + + + + - . |      |   | RMSE              |
          | + + + + + + + - . |      | X | Best RMSE (.0164) |
          | + + - - - - - - . |      | # | .0164 to .0173    |
          | + - + + - - - . . |      | + | .0173 to .0191    |
          | - - - - - - - . . |      | - | .0191 to .0208    |
        0 | - - - - - - - . . |      | . | .0208 to .0217    |
          +-------------------+      +-----------------------+
           0               2    Lambda unit

 Adapting search grid:
    Lambda unit  0 -------------[============]------- 2   [.8 ; 1.6]
    Lambda time  0 -----------[==========|============================================] 12   [2 ; 12]
    Lambda nn    0 []-------------------------------- 1   [0 ; .01] (and inf)

  Phase 2: Evaluating 10x10x11 (unit, time, nuclear norm) grid.

  First replication 1:53, roughly 3:08:20 expected in total
    5%  (5/100)   5:09 elapsed   about 1:37:51 remaining
   10%  (10/100)   5:54 elapsed   about 53:06 remaining
   15%  (15/100)   11:01 elapsed   about 1:02:25 remaining
   20%  (20/100)   11:45 elapsed   about 47:00 remaining
   25%  (25/100)   17:39 elapsed   about 52:57 remaining
   30%  (30/100)   18:30 elapsed   about 43:10 remaining
   35%  (35/100)   24:34 elapsed   about 45:37 remaining
   40%  (40/100)   25:23 elapsed   about 38:04 remaining
   45%  (45/100)   31:28 elapsed   about 38:27 remaining
   50%  (50/100)   32:20 elapsed   about 32:19 remaining
   55%  (55/100)   38:37 elapsed   about 31:35 remaining
   60%  (60/100)   39:28 elapsed   about 26:18 remaining
   65%  (65/100)   45:37 elapsed   about 24:33 remaining
   70%  (70/100)   46:29 elapsed   about 19:55 remaining
   75%  (75/100)   52:43 elapsed   about 17:34 remaining
   80%  (80/100)   53:37 elapsed   about 13:24 remaining
   85%  (85/100)   59:52 elapsed   about 10:33 remaining
   90%  (90/100)   1:00:49 elapsed   about 6:45 remaining
   95%  (95/100)   1:06:50 elapsed   about 3:31 remaining
  100%  (100/100)   1:07:35 elapsed
  Phase 2 chose lambdas (1.6, 3.11, .00667), on the edge of the searched range
    Lambda unit  .8 [===============================*] 1.6
    Lambda time  2 [======*=====|=====================================================] 12
    Lambda nn    0 [================|====*==========] .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
       12 | + + + + + + + + + + |
          | + + + + + + + + + + |
          | + + + + + + + + + + |
          | + + # + + # + + + # |      +- Legend --------------+
          | + + + + + + + + + + |      |   | RMSE              |
          | + + + + + + + + + + |      | X | Best RMSE (.0162) |
          | # # # # # # # # # # |      | # | .0162 to .0171    |
          | # # # # # # # # # # |      | + | .0171 to .019     |
          | # # # # # # # # # X |      | - | .019 to .0208     |
        2 | + + + + + + + + # # |      | . | .0208 to .0217    |
          +---------------------+      +-----------------------+
           .8               1.6    Lambda unit

 Expanding search grid: Lambda unit sits on the edge of the range.

  Phase 3: Evaluating 10x10x11 (unit, time, nuclear norm) grid.

   55%  (55/100)   6:34 elapsed   about 5:22 remaining
   60%  (60/100)   7:18 elapsed   about 4:52 remaining
   65%  (65/100)   13:53 elapsed   about 7:28 remaining
   70%  (70/100)   14:47 elapsed   about 6:20 remaining
   75%  (75/100)   21:23 elapsed   about 7:07 remaining
   80%  (80/100)   22:17 elapsed   about 5:34 remaining
   85%  (85/100)   28:51 elapsed   about 5:05 remaining
   90%  (90/100)   29:37 elapsed   about 3:17 remaining
   95%  (95/100)   36:01 elapsed   about 1:53 remaining
  100%  (100/100)   36:52 elapsed
  Phase 3 chose lambdas (1.87, 4.22, .00111)
    Lambda unit  .8 [===============|=====*==========] 2.4
    Lambda time  2 [============|=*===================================================] 12
    Lambda nn    0 [===*============|===============] .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
       12 | + + + + + + + + + + |
          | + + + + + + + + + + |
          | + + + + + + + + + + |
          | + + + + + # + + + + |      +- Legend --------------+
          | + + + + + + + + + + |      |   | RMSE              |
          | + + + + + + + + + + |      | X | Best RMSE (.0161) |
          | # # # # # # # # # + |      | # | .0161 to .017     |
          | # # # # # # X # # # |      | + | .017 to .0189     |
          | # # # # # # # # # # |      | - | .0189 to .0208    |
        2 | + + + + + + + + + + |      | . | .0208 to .0217    |
          +---------------------+      +-----------------------+
           .8               2.4    Lambda unit

 Adapting search grid:
    Lambda unit  .8 --------------[===============]--- 2.4   [1.47 ; 2.27]
    Lambda time  2 [===============================]----------------------------------- 12   [2 ; 6.72]
    Lambda nn    0 [===========]--------------------- .01   [0 ; .00361] (and inf)

  Evaluating 81 (unit x time) pairs in range, 10 nn values each
  First replication 2:22, roughly 3:11:42 expected in total
    5%  (5/81)   9:33 elapsed   about 2:25:09 remaining
   10%  (9/81)   12:31 elapsed   about 1:40:08 remaining
   15%  (13/81)   20:28 elapsed   about 1:47:03 remaining
   20%  (17/81)   24:29 elapsed   about 1:32:10 remaining
   25%  (21/81)   31:33 elapsed   about 1:30:08 remaining
   30%  (25/81)   36:28 elapsed   about 1:21:41 remaining
   35%  (29/81)   41:58 elapsed   about 1:15:15 remaining
   40%  (33/81)   48:15 elapsed   about 1:10:10 remaining
   45%  (37/81)   52:23 elapsed   about 1:02:17 remaining
   50%  (41/81)   59:27 elapsed   about 58:00 remaining
   55%  (45/81)   1:02:16 elapsed   about 49:48 remaining
   60%  (49/81)   1:10:38 elapsed   about 46:07 remaining
   65%  (53/81)   1:14:36 elapsed   about 39:24 remaining
   70%  (57/81)   1:21:56 elapsed   about 34:29 remaining
   75%  (61/81)   1:26:29 elapsed   about 28:21 remaining
   80%  (65/81)   1:34:58 elapsed   about 23:22 remaining
   85%  (69/81)   1:46:35 elapsed   about 18:32 remaining
   90%  (73/81)   1:53:07 elapsed   about 12:23 remaining
   95%  (77/81)   2:06:08 elapsed   about 6:33 remaining
  100%  (81/81)   2:10:30 elapsed

  Best after zoom (2.27, 4.95, .000451), 6:51:35 total
    Lambda unit  .8 --------------[==============*]--- 2.4
    Lambda time  2 [===================*===========]----------------------------------- 12
    Lambda nn    0 [*==========]--------------------- .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched in the adapted range.

Lambda time
     6.72 | . . . . . . . . . |
          | + + + - - + - - - |
          | + + + + + # + # + |      +- Legend -------------+
          | + + + + # # # # X |      |   | RMSE             |
          | + + + # # + # + + |      | X | Best RMSE (.016) |
          | + # # # # # # # # |      | # | .016 to .0162    |
          | + + # + + + + + # |      | + | .0162 to .0167   |
          | . . - - - - - - - |      | - | .0167 to .0171   |
        2 | . . . . . . - - - |      | . | .0171 to .0173   |
          +-------------------+      +----------------------+
          1.47            2.27    Lambda unit


Selected lambda: unit, time and nn = [2.266667 ; 4.951389 ; .0004514]

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.02145
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |           1
-------------+--------------------------------------------------
 lambda_unit |      2.2667
 lambda_time |      4.9514
   lambda_nn |    .0004514
             |  (selected by resample CV)
----------------------------------------------------------------
```

</details>

```
----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.02145
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |           1
-------------+--------------------------------------------------
 lambda_unit |      2.2667
 lambda_time |      4.9514
   lambda_nn |    .0004514
             |  (selected by resample CV)
----------------------------------------------------------------
```

Adaptive lands on lambdas close to what LOOCV picked (large `lambda_unit` and `lambda_time`, small `lambda_nn`), and reaches a lower CV RMSE than either cycle run (.016 against .0164 on the starting grid). It costs about seven hours here against twenty minutes for cycle, so it is worth reaching for when you suspect the cycle search has stopped at a local optimum, or when the chosen lambdas sit on the edge of the default grids, as they do under LOOCV above.

### Block adoption
 
Under block adoption, units are treated simultaneously in the final 18 periods. Under `group(time)`, which is the default, this is one pooled spell of 270 cells (15 units × 18 periods):

```s
trop y unit time w_block, cv(resample, seed(1))
```

which returns 

```
Cross-validating lambdas using resample with cycle search, and 200 trials (seed 1).
To reduce resample computational time, reduce no of trials or set lambdas.
  cycle 1 of up to 50:  lambda_unit -> 0   lambda_time -> 0   lambda_nn -> .1   5:58 elapsed
  cycle 2 of up to 50:  lambda_unit -> 0   lambda_time -> 0   lambda_nn -> .1   18:58 elapsed
  converged after 2 cycle(s), 18:58 total

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.01046
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          15
-------------+--------------------------------------------------
 lambda_unit |      0.0000
 lambda_time |      0.0000
   lambda_nn |          .1
             |  (selected by resample CV)
----------------------------------------------------------------
```

Under `group(cell)`, the same pattern produces 270 separate effects:
 
```s
trop y unit time w_block, group(cell) cv(loocv, cells(200) seed(1))
```
 
```
Cross-validating lambdas using loocv with cycle search, and 200 samples of 5058 total control cells.
To reduce loocv computational time, reduce number of cells or set lambdas.
  marginal:  lambda_time -> 4   lambda_unit -> 2   lambda_nn -> .25   1:37 elapsed
  cycle 1 of up to 50:  lambda_unit -> 2   lambda_time -> 4   lambda_nn -> .025   5:33 elapsed
  cycle 2 of up to 50:  lambda_unit -> 2   lambda_time -> 4   lambda_nn -> .025   16:28 elapsed
  converged after 2 cycle(s), 16:28 total
Computing 270 group effects.
  100%  (270/270)   0:04 elapsed

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |    -0.23351
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          15
-------------+--------------------------------------------------
 lambda_unit |      2.0000
 lambda_time |      4.0000
   lambda_nn |        .025
             |  (selected by loocv CV)
----------------------------------------------------------------
(270 per-cell effects: e(group_grid) [unit x time], e(group_tau), e(group_info))
```
 
The treatment effect are a lot higher using the parameters chosen by LOOCV. We'll study heterogenous treatment effects more closely by using `e(group_grid)`:
 
```s
matrix list e(group_grid), format(%9.4f)
```
 
```
e(group_grid)[15,18]
          t31      t32      t33      t34      t35      t36      t37      t38      t39      t40      t41      t42      t43      t44
 u97  -0.0019   0.0033  -0.0026  -0.0496  -0.1380  -0.9494  -1.0950  -1.1172  -1.1144  -1.1134  -1.1222  -1.0684  -0.9989  -0.9882
 u98   0.0233   0.0063  -0.0426  -0.2155  -0.2160  -0.8222  -0.9198  -0.9473  -0.9809  -0.9648  -0.9726  -0.9572  -0.9014  -0.9146
 u99   0.0715   0.1253   0.1810   0.2353   0.2291  -0.0618  -0.1138  -0.1597  -0.2089  -0.1765  -0.1540  -0.1198  -0.0469  -0.0502
u100  -0.0079  -0.0014  -0.0321  -0.0560  -0.0522   0.0158   0.0359   0.0387   0.0603   0.1105   0.1458   0.1856   0.2672   0.2926
u101   0.0316   0.0408   0.0725   0.0636   0.0188  -0.2378  -0.2742  -0.2865  -0.2820  -0.2435  -0.2214  -0.1813  -0.1264  -0.1366
(output truncated)
```
 
Under block adoption, LOOCV with `group(cell)` performs poorly relative to resample with `group(time)`. The heterogenous treatment effects above shows that units further from the onset of treatment have larger treatment effects. This is because LOOCV estimates placebo effects using control units which are surrounded by adjacent donor units. CV minimises the criterion under this design, and chooses lambdas which are large (`lambda_time=4`). This heavily weights nearby units. However, treated units in block actually sit far from control units (i.e. up to 19 periods). Resample CV under `group(time)` suits block adoption, because the CV uses the actual treatment pattern on a subset of never-treated control units from the panel. As a result, we recommend using resample under block adoption with many treated periods. 
 
### Staggered adoption
 
Next, we study staggered adoption with three cohorts of five units. Under `group(time)`, the `detail` option prints a summary which includes which units are included in each cohort:
 
```s
trop y unit time w_stag, cv(resample, seed(1)) detail
```
 
```
Cross-validating lambdas using resample with cycle search, and 200 trials (seed 1).
To reduce resample computational time, reduce no of trials or set lambdas.
  cycle 1 of up to 50:  lambda_unit -> 0   lambda_time -> 0   lambda_nn -> .005   21:37 elapsed
  cycle 2 of up to 50:  lambda_unit -> 0   lambda_time -> 0   lambda_nn -> .005   2:42:31 elapsed
  converged after 2 cycle(s), 2:42:31 total

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |    -0.01011
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          15
-------------+--------------------------------------------------
 lambda_unit |      0.0000
 lambda_time |      0.0000
   lambda_nn |        .005
             |  (selected by resample CV)
----------------------------------------------------------------
(3 per-spell effects in e(group_tau), e(group_weight), e(group_info), e(group_units))

time periods  start   end  units  cells        tau   cohort's units
      t21_48     21    48      5    140    -0.0453   97 98 99 100 101
      t31_48     31    48      5     90     0.0595   102 103 104 105 106
      t41_48     41    48      5     40    -0.0435   107 108 109 110 111
```
 
The per-cohort effects are also returned in `e(group_tau)` if you have not used the detail option (but the indexed cohort units are not):
 
```s
matrix list e(group_tau)
```
 
```
e(group_tau)[1,3]
        t21_48      t31_48      t41_48
r1  -.04527924   .05947124  -.04354185
```
### General assignment
 
Under general assignment, treatment can switch on and off. Under `group(time)`, treated periods with the same start and end are pooled across units and single units remain their own estimand. Here the treated units generate 51 individual estimands:
 
```s
trop y unit time w_gen, lambda_unit(0.3) lambda_time(0.5) lambda_nn(0.025)
```
 
```
Computing 51 group effects.
  100%  (51/51)   0:22 elapsed

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.00128
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          10
-------------+--------------------------------------------------
 lambda_unit |      0.3000
 lambda_time |      0.5000
   lambda_nn |        .025
----------------------------------------------------------------
(51 per-spell effects in e(group_tau), e(group_weight), e(group_info), e(group_units))
```
 
`e(group_info)` records each estimands start and end time, as well as units:
 
```s
matrix list e(group_info)
```
 
```
e(group_info)[51,4]
          start      end  n_units  n_cells
t10_12       10       12        1        3
t18_18       18       18        2        2
t22_22       22       22        1        1
t28_28       28       28        2        2
t33_35       33       35        1        3
t39_40       39       40        1        2
t45_45       45       45        1        1
  t1_1        1        1        2        2
(output truncated)
```

### Inference

A stratified block bootstrap gives a standard error and a percentile confidence
interval. Tune once and bootstrap at the selected regularizers (rather than
re-tuning inside every replication):

```s
trop y unit time w, lambda_unit(0) lambda_time(1) lambda_nn(0.1) vce(bootstrap)
```

which returns

```
To reduce computational time, reduce reps() or use vce(jackknife).
Bootstrap inference using 200 bootstrap replications.
  100%  (200/200)   0:13 elapsed

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.02963
   Std. err. |     0.02087
      95% CI |   -0.01074    0.06781
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          29
-------------+--------------------------------------------------
 lambda_unit |      0.0000
 lambda_time |      1.0000
   lambda_nn |          .1
----------------------------------------------------------------
```
Next, we also provide an option for jackknife inference which also shows a null effect.

```s
trop y unit time w, lambda_unit(0) lambda_time(1) lambda_nn(0.1) vce(jackknife)
```

```
Jackknife inference using 111 leave-one-out replications.
  100%  (111/111)   0:09 elapsed

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.02963
   Std. err. |     0.02242
      95% CI |   -0.01431    0.07357
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          29
-------------+--------------------------------------------------
 lambda_unit |      0.0000
 lambda_time |      1.0000
   lambda_nn |          .1
----------------------------------------------------------------
```

## Computational costs
The examples in this readme were run on the following computational setup:
1.	Operating System:    Windows 11 Enterprise (Version 10.0.22621 Build 22621)
2.	CPU:   13th Gen Intel(R) Core(TM) i7-1355U, 1700 Mhz, 10 Core(s), 12 Logical Processor(s)
3.	RAM:   32 GB
4.	Software version: Stata/SE 18.0 for Windows (64-bit x86-64)

In general, we advise:
1. Testing your data using TROP with set lambdas for fast computations to begin with.
2. Running CV on the main specification (resample best suited for designs with long, uninterrupted treatment periods).
3. Running block bootstrap and jackknife inference using the lambdas from cross-validation.

## References
Abadie, A., Diamond, A., & Hainmueller, J. (2010). [Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program](https://doi.org/10.1198/jasa.2009.ap08746). *Journal of the American Statistical Association*, 105(490), 493–505.

Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021). [Synthetic difference-in-differences](https://doi.org/10.1257/aer.20190159). *American Economic Review*, 111(12), 4088–4118.

Athey, S., Imbens, G., Qu, Z., & Viviano, D. (2025). [Triply robust panel estimators](https://arxiv.org/pdf/2508.21536). arXiv preprint arXiv:2508.21536.
