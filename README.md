# TROP: Triply Robust Panel Estimators

This Stata package implements Triply Robust Panel Estimators (TROP) following [Athey et al. (2025)](#references). TROP is a panel-data estimator for treatment effects. It combines a flexible model for potential outcomes, with unit weights which upweight units similar to the treated units, and time weights which upweight the time periods closest to treatment.

The latest version of `trop` can be installed directly from GitHub in Stata:

```s
net install trop, from("https://raw.githubusercontent.com/justinwaddy/TROP/main/") replace
```

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

No unit is actually treated in this null effect study, so the true effect is zero and the estimate here is small. Using set lambdas are good for a first run, but they were chosen arbitrarily. Let's use cross-validation to choose the lambdas using a data-driven approach instead. We will begin by using resampling cross-validation:
```s
trop y unit time w_single, cv(resample, seed(1))
```
 
```
Cross-validating lambdas using resample with cycle search, and 200 trials (seed 1).
To reduce resample computational time, reduce no of trials or set lambdas.
  cycle 1 of up to 50:  lambda_unit -> .2   lambda_time -> 2   lambda_nn -> .05   4:38 elapsed
  cycle 2 of up to 50:  lambda_unit -> .8   lambda_time -> 2   lambda_nn -> .05   12:24 elapsed
  cycle 3 of up to 50:  lambda_unit -> .8   lambda_time -> 2   lambda_nn -> .05   20:05 elapsed
  converged after 3 cycle(s), 20:05 total

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
 
This takes about twenty minutes and moves the estimate slightly closer to zero. The cycle search converged after three passes at `lambda_unit = 0.8` and `lambda_time = 2`, both interior to the default grids. We used cycle but it may be that our lambda grid which we searched was too sparse. The `adaptive` search instead sweeps a coarse grid, widens it whenever the winner lands on an edge, and zooms in on the best point. With `verbose` it prints an RMSE heat map after each phase which you can see by expanding the output below:

```s
trop y unit time w_single, cv(resample adaptive, trials(200)) verbose
```

<details>
<summary><b>Full adaptive cross-validation log</b> (4 phases, 11 hours 37 minutes total)</summary>

```
Cross-validating lambdas using resample with adaptive search, and 200 trials (seed 0).
To reduce resample computational time, reduce no of trials or set lambdas.
Adaptive CV: phase 1 sweeps the 9x11x9 starting grid of lambdas;
later phases refine around the winner; points() and expansions() control their cost.
 
  Phase 1: Sweeping the 9x11x9 (unit, time, nuclear norm) starting grid.

  First replication 0:35, roughly 57:45 expected in total
    5%  (5/99)   9:21 elapsed   about 2:55:46 remaining
   10%  (10/99)   22:48 elapsed   about 3:22:55 remaining
   15%  (15/99)   31:43 elapsed   about 2:57:36 remaining
   20%  (20/99)   49:04 elapsed   about 3:13:48 remaining
   25%  (25/99)   55:15 elapsed   about 2:43:32 remaining
   30%  (30/99)   1:03:48 elapsed   about 2:26:44 remaining
   35%  (35/99)   1:07:33 elapsed   about 2:03:31 remaining
   40%  (40/99)   1:14:56 elapsed   about 1:50:31 remaining
   45%  (45/99)   1:20:26 elapsed   about 1:36:31 remaining
   50%  (50/99)   1:27:56 elapsed   about 1:26:10 remaining
   55%  (55/99)   1:34:50 elapsed   about 1:15:52 remaining
   60%  (60/99)   1:49:51 elapsed   about 1:11:24 remaining
   65%  (65/99)   2:03:08 elapsed   about 1:04:24 remaining
   70%  (70/99)   2:11:58 elapsed   about 54:40 remaining
   75%  (75/99)   2:21:56 elapsed   about 45:25 remaining
   80%  (80/99)   2:29:02 elapsed   about 35:23 remaining
   85%  (85/99)   2:37:02 elapsed   about 25:51 remaining
   90%  (90/99)   2:45:32 elapsed   about 16:33 remaining
   95%  (95/99)   2:57:34 elapsed   about 7:28 remaining
  100%  (99/99)   2:59:26 elapsed
  Phase 1 chose lambdas (1.2, 4, .005), on the edge of the swept grid
    Lambda unit  0 .1 .2 .3 .5 .8 [1.2] 1.6 2
    Lambda time  0 .025 .05 .1 .2 .35 .5 .75 1 2 [4]
    Lambda nn    [.005] .01 .025 .05 .1 .25 .5 1 inf

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
        4 | + # # # # # X # # |
          | # # # # # # # # # |
          | + + + + + + + + # |
          | + + + + + + + # # |
          | + + + + + + + + - |      +- Legend --------------+
          | + + + + + + + - - |      |   | RMSE              |
          | + + + + + + + - . |      | X | Best RMSE (.0165) |
          | + + - - - - - - . |      | # | .0165 to .0174    |
          | + + + + - - - . . |      | + | .0174 to .0191    |
          | - - - - - - - . . |      | - | .0191 to .0208    |
        0 | - - - - - - - . . |      | . | .0208 to .0217    |
          +-------------------+      +-----------------------+
           0               2    Lambda unit

 Adapting search grid:
    Lambda unit  0 -------------[============]------- 2   [.8 ; 1.6]
    Lambda time  0 -----------[==========|============================================] 12   [2 ; 12]
    Lambda nn    0 []-------------------------------- 1   [0 ; .01] (and inf)

  Phase 2: Evaluating 10x10x11 (unit, time, nuclear norm) grid.

  First replication 1:54, roughly 3:10:00 expected in total
    5%  (5/100)   5:06 elapsed   about 1:36:54 remaining
   10%  (10/100)   6:36 elapsed   about 59:24 remaining
   15%  (15/100)   20:54 elapsed   about 1:58:25 remaining
   20%  (20/100)   23:43 elapsed   about 1:34:52 remaining
   25%  (25/100)   40:11 elapsed   about 2:00:33 remaining
   30%  (30/100)   43:10 elapsed   about 1:40:43 remaining
   35%  (35/100)   58:25 elapsed   about 1:48:29 remaining
   40%  (40/100)   1:00:56 elapsed   about 1:31:24 remaining
   45%  (45/100)   1:15:02 elapsed   about 1:31:42 remaining
   50%  (50/100)   1:17:32 elapsed   about 1:17:32 remaining
   55%  (55/100)   1:32:13 elapsed   about 1:15:27 remaining
   60%  (60/100)   1:34:40 elapsed   about 1:03:06 remaining
   65%  (65/100)   1:49:07 elapsed   about 58:45 remaining
   70%  (70/100)   1:51:43 elapsed   about 47:52 remaining
   75%  (75/100)   2:05:47 elapsed   about 41:55 remaining
   80%  (80/100)   2:08:15 elapsed   about 32:03 remaining
   85%  (85/100)   2:22:26 elapsed   about 25:08 remaining
   90%  (90/100)   2:24:49 elapsed   about 16:05 remaining
   95%  (95/100)   2:38:59 elapsed   about 8:22 remaining
  100%  (100/100)   2:41:17 elapsed
  Phase 2 chose lambdas (1.6, 4.22, .00222), on the edge of the searched range
    Lambda unit  .8 [===============================*] 1.6
    Lambda time  2 [============|=*===================================================] 12
    Lambda nn    0 [======*=========|===============] .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
       12 | + + + + + + + + + + |
          | + + + + + + + + + + |
          | + + + + + + + + + + |
          | + + + # # + + + + + |      +- Legend --------------+
          | + + + + + + + + + + |      |   | RMSE              |
          | + + + + + + + + + + |      | X | Best RMSE (.0162) |
          | # # # # # # # # # # |      | # | .0162 to .0171    |
          | # # # # # # # # # X |      | + | .0171 to .019     |
          | # # # # # # # # # # |      | - | .019 to .0208     |
        2 | + + + + + + + + # # |      | . | .0208 to .0217    |
          +---------------------+      +-----------------------+
           .8               1.6    Lambda unit

 Expanding search grid: Lambda unit sits on the edge of the range.

  Phase 3: Evaluating 10x10x11 (unit, time, nuclear norm) grid.

   55%  (55/100)   14:06 elapsed   about 11:32 remaining
   60%  (60/100)   16:28 elapsed   about 10:58 remaining
   65%  (65/100)   30:51 elapsed   about 16:36 remaining
   70%  (70/100)   33:17 elapsed   about 14:15 remaining
   75%  (75/100)   47:13 elapsed   about 15:44 remaining
   80%  (80/100)   49:39 elapsed   about 12:24 remaining
   85%  (85/100)   1:03:25 elapsed   about 11:11 remaining
   90%  (90/100)   1:05:52 elapsed   about 7:19 remaining
   95%  (95/100)   1:19:21 elapsed   about 4:10 remaining
  100%  (100/100)   1:21:54 elapsed
  Phase 3 chose lambdas (2.04, 4.22, .00111)
    Lambda unit  .8 [===============|=========*======] 2.4
    Lambda time  2 [============|=*===================================================] 12
    Lambda nn    0 [===*============|===============] .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
       12 | + + + + + + + + + + |
          | + + + + + + + + + + |
          | + + + + + + + + + + |
          | + + # + + + + + + + |      +- Legend --------------+
          | + + + + + # + + + + |      |   | RMSE              |
          | + + + + + + + + + + |      | X | Best RMSE (.0162) |
          | # # # # # # # # # # |      | # | .0162 to .0171    |
          | # # # # # # # X # # |      | + | .0171 to .0189    |
          | # # # # # # # # # # |      | - | .0189 to .0208    |
        2 | + + + + + + + # # # |      | . | .0208 to .0217    |
          +---------------------+      +-----------------------+
           .8               2.4    Lambda unit

 Adapting search grid:
    Lambda unit  .8 -----------------[===============] 2.4   [1.64 ; 2.4]
    Lambda time  2 [===============================]----------------------------------- 12   [2 ; 6.72]
    Lambda nn    0 [===========]--------------------- .01   [0 ; .00361] (and inf)

  Evaluating 81 (unit x time) pairs in range, 10 nn values each
  First replication 5:56, roughly 8:00:36 expected in total
    5%  (5/81)   23:08 elapsed   about 5:51:37 remaining
   10%  (9/81)   30:00 elapsed   about 4:00:00 remaining
   15%  (13/81)   50:18 elapsed   about 4:23:06 remaining
   20%  (17/81)   59:53 elapsed   about 3:45:26 remaining
   25%  (21/81)   1:18:02 elapsed   about 3:42:57 remaining
   30%  (25/81)   1:30:11 elapsed   about 3:22:00 remaining
   35%  (29/81)   1:45:32 elapsed   about 3:09:13 remaining
   40%  (33/81)   2:00:11 elapsed   about 2:54:48 remaining
   45%  (37/81)   2:10:53 elapsed   about 2:35:38 remaining
   50%  (41/81)   2:27:52 elapsed   about 2:24:15 remaining
   55%  (45/81)   2:34:17 elapsed   about 2:03:25 remaining
   60%  (49/81)   2:53:12 elapsed   about 1:53:06 remaining
   65%  (53/81)   3:02:29 elapsed   about 1:36:24 remaining
   70%  (57/81)   3:19:41 elapsed   about 1:24:04 remaining
   75%  (61/81)   3:30:41 elapsed   about 1:09:04 remaining
   80%  (65/81)   3:44:52 elapsed   about 55:21 remaining
   85%  (69/81)   3:59:08 elapsed   about 41:35 remaining
   90%  (73/81)   4:10:01 elapsed   about 27:23 remaining
   95%  (77/81)   4:27:55 elapsed   about 13:55 remaining
  100%  (81/81)   4:34:51 elapsed

  Best after zoom (2.12, 4.36, .000903), 11:37:28 total
    Lambda unit  .8 -----------------[=========*=====] 2.4
    Lambda time  2 [===============*===============]----------------------------------- 12
    Lambda nn    0 [==*========]--------------------- .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched in the adapted range.

Lambda time
     6.72 | . . . . . . . . . |
          | + + + - - - - - - |
          | + + + + + + + + + |      +- Legend -------------+
          | + + + # # # # # + |      |   | RMSE             |
          | # # # + # X # # + |      | X | Best RMSE (.016) |
          | # # # # # # # # # |      | # | .016 to .0163    |
          | + # + + + + + # + |      | + | .0163 to .0167   |
          | . . - - - - - - - |      | - | .0167 to .0171   |
        2 | . . . . . . . - - |      | . | .0171 to .0173   |
          +-------------------+      +----------------------+
          1.64            2.4    Lambda unit


Selected lambda: unit, time and nn = [2.116667 ; 4.361111 ; .0009028]

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.02242
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |           1
-------------+--------------------------------------------------
 lambda_unit |      2.1167
 lambda_time |      4.3611
   lambda_nn |    .0009028
             |  (selected by resample CV)
----------------------------------------------------------------
```

</details>

```
----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.02242
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |           1
-------------+--------------------------------------------------
 lambda_unit |      2.1167
 lambda_time |      4.3611
   lambda_nn |    .0009028
             |  (selected by resample CV)
----------------------------------------------------------------
```

Adaptive reaches a lower CV RMSE than the cycle. It costs about twelve hours here against twenty minutes for cycle, so it is worth reaching for when you suspect the cycle search has stopped at a local optimum or when the chosen lambdas sit on the edge of the default grids.

Next, let's try leave-one-out cross validation (LOOCV):
```s
trop y unit time w_single, group(cell) cv(loocv, cells(200) seed(1))
```
 
```
Cross-validating lambdas using loocv with cycle search, and 200 samples of 5327 total control cells.
To reduce loocv computational time, reduce number of cells or set lambdas.
  marginal:  lambda_time -> .2   lambda_unit -> 2   lambda_nn -> .005   0:40 elapsed
  cycle 1 of up to 50:  lambda_unit -> .2   lambda_time -> .2   lambda_nn -> .25   1:54 elapsed
  cycle 2 of up to 50:  lambda_unit -> .2   lambda_time -> .2   lambda_nn -> .25   2:59 elapsed
  converged after 2 cycle(s), 2:59 total

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.02725
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |           1
-------------+--------------------------------------------------
 lambda_unit |      0.2000
 lambda_time |      0.2000
   lambda_nn |         .25
             |  (selected by loocv CV)
----------------------------------------------------------------
```

LOOCV finished in 2:59 but scored only 200 of the 5,327 control cells, set by `cells(200)`. Full LOOCV on this panel is considerably slower than the resample run. LOOCV lands slightly further from zero (0.027 against 0.022) and chooses much weaker unit and time weights, `lambda_unit = 0.2` and `lambda_time = 0.2`, both interior to the default grids. The cycle search only moves one lambda at a time over a coarse grid, so let's try LOOCV using the adaptive approach instead:

```s
trop y unit time w_single, group(cell) cv(loocv adaptive, cells(200)) verbose
```

<details>
<summary><b>Full adaptive cross-validation log</b> (6 phases plus zoom, 2 hours 12 minutes total)</summary>

```
Cross-validating lambdas using loocv with adaptive search, and 200 samples of 5327 total control cells.
To reduce loocv computational time, reduce number of cells or set lambdas.
Adaptive CV: phase 1 sweeps the 9x11x9 starting grid of lambdas;
later phases refine around the winner; points() and expansions() control their cost.
 
  Phase 1: Sweeping the 9x11x9 (unit, time, nuclear norm) starting grid.

  First replication 0:11, roughly 18:09 expected in total
    5%  (5/99)   0:58 elapsed   about 18:10 remaining
   10%  (10/99)   2:32 elapsed   about 22:32 remaining
   15%  (15/99)   3:41 elapsed   about 20:37 remaining
   20%  (20/99)   5:01 elapsed   about 19:48 remaining
   25%  (25/99)   6:17 elapsed   about 18:35 remaining
   30%  (30/99)   7:26 elapsed   about 17:05 remaining
   35%  (35/99)   8:52 elapsed   about 16:12 remaining
   40%  (40/99)   9:51 elapsed   about 14:31 remaining
   45%  (45/99)   11:26 elapsed   about 13:43 remaining
   50%  (50/99)   12:23 elapsed   about 12:08 remaining
   55%  (55/99)   14:03 elapsed   about 11:14 remaining
   60%  (60/99)   15:03 elapsed   about 9:46 remaining
   65%  (65/99)   16:39 elapsed   about 8:42 remaining
   70%  (70/99)   17:43 elapsed   about 7:20 remaining
   75%  (75/99)   18:59 elapsed   about 6:04 remaining
   80%  (80/99)   20:16 elapsed   about 4:48 remaining
   85%  (85/99)   21:22 elapsed   about 3:31 remaining
   90%  (90/99)   22:37 elapsed   about 2:15 remaining
   95%  (95/99)   23:36 elapsed   about 0:59 remaining
  100%  (99/99)   24:42 elapsed
  Phase 1 chose lambdas (0, 0, .005), on the edge of the swept grid
    Lambda unit  [0] .1 .2 .3 .5 .8 1.2 1.6 2
    Lambda time  [0] .025 .05 .1 .2 .35 .5 .75 1 2 4
    Lambda nn    [.005] .01 .025 .05 .1 .25 .5 1 inf

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
        4 | . - - - . . . . . |
          | - - - - - - - - - |
          | + + - + + - + + + |
          | - + + + + + + + + |
          | # # # # # # # # # |      +- Legend -----------------+
          | # # # # # # # # # |      |   | RMSE                 |
          | # # # # # # # # # |      | X | Best RMSE (3.86e-12) |
          | # # # # # # # # # |      | # | 3.86e-12 to 1.23e-06 |
          | # # # # # # # # # |      | + | 1.23e-06 to 3.70e-06 |
          | # # # # # # # # # |      | - | 3.70e-06 to 6.17e-06 |
        0 | X # # # # # # # # |      | . | 6.17e-06 to 7.41e-06 |
          +-------------------+      +--------------------------+
           0               2    Lambda unit

 Adapting search grid:
    Lambda unit  0 [=]------------------------------- 2   [0 ; .1]
    Lambda time  0 []-------------------------------- 4   [0 ; .025]
    Lambda nn    0 []-------------------------------- 1   [0 ; .01] (and inf)

  Phase 2: Evaluating 10x10x11 (unit, time, nuclear norm) grid.

  First replication 0:14, roughly 23:20 expected in total
    5%  (5/100)   1:09 elapsed   about 21:51 remaining
   10%  (10/100)   2:18 elapsed   about 20:42 remaining
   15%  (15/100)   3:30 elapsed   about 19:50 remaining
   20%  (20/100)   4:41 elapsed   about 18:44 remaining
   25%  (25/100)   5:51 elapsed   about 17:33 remaining
   30%  (30/100)   7:02 elapsed   about 16:24 remaining
   35%  (35/100)   8:11 elapsed   about 15:11 remaining
   40%  (40/100)   9:24 elapsed   about 14:06 remaining
   45%  (45/100)   10:34 elapsed   about 12:54 remaining
   50%  (50/100)   11:44 elapsed   about 11:44 remaining
   55%  (55/100)   12:54 elapsed   about 10:33 remaining
   60%  (60/100)   14:04 elapsed   about 9:22 remaining
   65%  (65/100)   15:15 elapsed   about 8:12 remaining
   70%  (70/100)   16:25 elapsed   about 7:02 remaining
   75%  (75/100)   17:35 elapsed   about 5:51 remaining
   80%  (80/100)   18:48 elapsed   about 4:42 remaining
   85%  (85/100)   19:59 elapsed   about 3:31 remaining
   90%  (90/100)   21:10 elapsed   about 2:21 remaining
   95%  (95/100)   22:21 elapsed   about 1:10 remaining
  100%  (100/100)   23:32 elapsed
  Phase 2 chose lambdas (.1, .025, inf), on the edge of the searched range
    Lambda unit  0 [===============================*] .1
    Lambda time  0 [===============================*] .025
    Lambda nn    0 [================|===============] .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
     .025 | # # # # # # # # # X |
          | # # # # # # # # # # |
          | # # # # # # # # # # |
          | # # # # # # # # # # |      +- Legend -----------------+
          | # # # # # # # # # # |      |   | RMSE                 |
          | # # # # # # # # # # |      | X | Best RMSE (1.95e-11) |
          | # # # # # # # # # # |      | # | 3.86e-12 to 1.23e-06 |
          | # # # # # # # # # # |      | + | 1.23e-06 to 3.70e-06 |
          | # # # # # # # # # # |      | - | 3.70e-06 to 6.17e-06 |
        0 | # # # # # # # # # # |      | . | 6.17e-06 to 7.41e-06 |
          +---------------------+      +--------------------------+
           0                 .1    Lambda unit

 Expanding search grid: Lambda unit and Lambda time sit on the edge of the range.

  Phase 3: Evaluating 10x10x11 (unit, time, nuclear norm) grid.

   10%  (10/100)   1:11 elapsed   about 10:39 remaining
   20%  (20/100)   2:24 elapsed   about 9:36 remaining
   30%  (30/100)   3:36 elapsed   about 8:24 remaining
   40%  (40/100)   4:48 elapsed   about 7:12 remaining
   50%  (50/100)   5:59 elapsed   about 5:59 remaining
   55%  (55/100)   7:11 elapsed   about 5:52 remaining
   60%  (60/100)   8:23 elapsed   about 5:35 remaining
   65%  (65/100)   9:34 elapsed   about 5:09 remaining
   70%  (70/100)   10:45 elapsed   about 4:36 remaining
   75%  (75/100)   11:56 elapsed   about 3:58 remaining
   80%  (80/100)   13:07 elapsed   about 3:16 remaining
   85%  (85/100)   14:19 elapsed   about 2:31 remaining
   90%  (90/100)   15:30 elapsed   about 1:43 remaining
   95%  (95/100)   16:40 elapsed   about 0:52 remaining
  100%  (100/100)   17:52 elapsed
  Phase 3 chose lambdas (.2, .05, inf), on the edge of the searched range
    Lambda unit  0 [================|==============*] .2
    Lambda time  0 [================|==============*] .05
    Lambda nn    0 [================|===============] .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
      .05 | # # # # # # # # # X |
          | # # # # # # # # # # |
          | # # # # # # # # # # |
          | # # # # # # # # # # |      +- Legend -----------------+
          | # # # # # # # # # # |      |   | RMSE                 |
          | # # # # # # # # # # |      | X | Best RMSE (1.57e-11) |
          | # # # # # # # # # # |      | # | 3.86e-12 to 1.23e-06 |
          | # # # # # # # # # # |      | + | 1.23e-06 to 3.70e-06 |
          | # # # # # # # # # # |      | - | 3.70e-06 to 6.17e-06 |
        0 | # # # # # # # # # # |      | . | 6.17e-06 to 7.41e-06 |
          +---------------------+      +--------------------------+
           0                 .2    Lambda unit

 Expanding search grid: Lambda unit and Lambda time sit on the edge of the range.

  Phase 4: Evaluating 10x10x11 (unit, time, nuclear norm) grid.

   10%  (10/100)   1:11 elapsed   about 10:39 remaining
   20%  (20/100)   2:22 elapsed   about 9:28 remaining
   30%  (30/100)   3:34 elapsed   about 8:19 remaining
   40%  (40/100)   4:45 elapsed   about 7:07 remaining
   50%  (50/100)   5:57 elapsed   about 5:57 remaining
   55%  (55/100)   7:10 elapsed   about 5:51 remaining
   60%  (60/100)   8:21 elapsed   about 5:34 remaining
   65%  (65/100)   9:33 elapsed   about 5:08 remaining
   70%  (70/100)   10:45 elapsed   about 4:36 remaining
   75%  (75/100)   11:57 elapsed   about 3:59 remaining
   80%  (80/100)   13:08 elapsed   about 3:17 remaining
   85%  (85/100)   14:22 elapsed   about 2:32 remaining
   90%  (90/100)   15:33 elapsed   about 1:43 remaining
   95%  (95/100)   16:44 elapsed   about 0:52 remaining
  100%  (100/100)   17:56 elapsed
  Phase 4 chose lambdas (.4, .1, inf), on the edge of the searched range
    Lambda unit  0 [================|==============*] .4
    Lambda time  0 [================|==============*] .1
    Lambda nn    0 [================|===============] .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
       .1 | # # # # # # # # # X |
          | # # # # # # # # # # |
          | # # # # # # # # # # |
          | # # # # # # # # # # |      +- Legend -----------------+
          | # # # # # # # # # # |      |   | RMSE                 |
          | # # # # # # # # # # |      | X | Best RMSE (1.08e-11) |
          | # # # # # # # # # # |      | # | 3.86e-12 to 1.23e-06 |
          | # # # # # # # # # # |      | + | 1.23e-06 to 3.70e-06 |
          | # # # # # # # # # # |      | - | 3.70e-06 to 6.17e-06 |
        0 | # # # # # # # # # # |      | . | 6.17e-06 to 7.41e-06 |
          +---------------------+      +--------------------------+
           0                 .4    Lambda unit

 Expanding search grid: Lambda unit and Lambda time sit on the edge of the range.

  Phase 5: Evaluating 10x10x11 (unit, time, nuclear norm) grid.

   10%  (10/100)   1:11 elapsed   about 10:39 remaining
   20%  (20/100)   2:20 elapsed   about 9:20 remaining
   30%  (30/100)   3:31 elapsed   about 8:12 remaining
   40%  (40/100)   4:45 elapsed   about 7:07 remaining
   50%  (50/100)   5:55 elapsed   about 5:55 remaining
   55%  (55/100)   7:08 elapsed   about 5:50 remaining
   60%  (60/100)   8:24 elapsed   about 5:36 remaining
   65%  (65/100)   9:37 elapsed   about 5:10 remaining
   70%  (70/100)   10:49 elapsed   about 4:38 remaining
   75%  (75/100)   12:01 elapsed   about 4:00 remaining
   80%  (80/100)   13:12 elapsed   about 3:18 remaining
   85%  (85/100)   14:24 elapsed   about 2:32 remaining
   90%  (90/100)   15:35 elapsed   about 1:43 remaining
   95%  (95/100)   16:48 elapsed   about 0:53 remaining
  100%  (100/100)   17:59 elapsed
  Phase 5 chose lambdas (.0889, .2, inf), on the edge of the searched range
    Lambda unit  0 [===*============|===============] .8
    Lambda time  0 [================|==============*] .2
    Lambda nn    0 [================|===============] .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
       .2 | # X # # # # # # # # |
          | # # # # # # # # # # |
          | # # # # # # # # # # |
          | # # # # # # # # # # |      +- Legend -----------------+
          | # # # # # # # # # # |      |   | RMSE                 |
          | # # # # # # # # # # |      | X | Best RMSE (7.12e-12) |
          | # # # # # # # # # # |      | # | 3.86e-12 to 1.23e-06 |
          | # # # # # # # # # # |      | + | 1.23e-06 to 3.70e-06 |
          | # # # # # # # # # # |      | - | 3.70e-06 to 6.17e-06 |
        0 | # # # # # # # # # # |      | . | 6.17e-06 to 7.41e-06 |
          +---------------------+      +--------------------------+
           0                 .8    Lambda unit

 Expanding search grid: Lambda time sits on the edge of the range.

  Phase 6: Evaluating 10x10x11 (unit, time, nuclear norm) grid.

   10%  (10/100)   1:12 elapsed   about 10:48 remaining
   20%  (20/100)   2:23 elapsed   about 9:32 remaining
   30%  (30/100)   3:35 elapsed   about 8:21 remaining
   40%  (40/100)   4:46 elapsed   about 7:09 remaining
   50%  (50/100)   5:58 elapsed   about 5:58 remaining
   60%  (60/100)   7:09 elapsed   about 4:46 remaining
   70%  (70/100)   8:21 elapsed   about 3:34 remaining
   80%  (80/100)   9:34 elapsed   about 2:23 remaining
   90%  (90/100)   10:47 elapsed   about 1:11 remaining
  100%  (100/100)   11:59 elapsed
  Phase 6 chose lambdas (.533, .178, inf)
    Lambda unit  0 [================|====*==========] .8
    Lambda time  0 [==============*=|===============] .4
    Lambda nn    0 [================|===============] .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched.

Lambda time
       .4 | # # # # # # # # # # |
          | # # # # # # # # # # |
          | # # # # # # # # # # |
          | # # # # # # # # # # |      +- Legend -----------------+
          | # # # # # # # # # # |      |   | RMSE                 |
          | # # # # # # X # # # |      | X | Best RMSE (7.38e-12) |
          | # # # # # # # # # # |      | # | 3.86e-12 to 1.23e-06 |
          | # # # # # # # # # # |      | + | 1.23e-06 to 3.70e-06 |
          | # # # # # # # # # # |      | - | 3.70e-06 to 6.17e-06 |
        0 | # # # # # # # # # # |      | . | 6.17e-06 to 7.41e-06 |
          +---------------------+      +--------------------------+
           0                 .8    Lambda unit

 Adapting search grid:
    Lambda unit  0 --------------[===============]--- .8   [.333 ; .733]
    Lambda time  0 ------[================]---------- .4   [.0778 ; .278]
    Lambda nn    0 [================================] .01   [0 ; .01] (and inf)

  Evaluating 81 (unit x time) pairs in range, 10 nn values each
  First replication 0:13, roughly 17:33 expected in total
    5%  (5/81)   1:08 elapsed   about 17:13 remaining
   10%  (9/81)   2:00 elapsed   about 16:00 remaining
   15%  (13/81)   2:54 elapsed   about 15:10 remaining
   20%  (17/81)   3:46 elapsed   about 14:10 remaining
   25%  (21/81)   4:39 elapsed   about 13:17 remaining
   30%  (25/81)   5:31 elapsed   about 12:21 remaining
   35%  (29/81)   6:25 elapsed   about 11:30 remaining
   40%  (33/81)   7:17 elapsed   about 10:35 remaining
   45%  (37/81)   8:10 elapsed   about 9:42 remaining
   50%  (41/81)   9:04 elapsed   about 8:50 remaining
   55%  (45/81)   9:57 elapsed   about 7:57 remaining
   60%  (49/81)   10:51 elapsed   about 7:05 remaining
   65%  (53/81)   11:42 elapsed   about 6:10 remaining
   70%  (57/81)   12:35 elapsed   about 5:17 remaining
   75%  (61/81)   13:32 elapsed   about 4:26 remaining
   80%  (65/81)   14:30 elapsed   about 3:34 remaining
   85%  (69/81)   15:24 elapsed   about 2:40 remaining
   90%  (73/81)   16:16 elapsed   about 1:46 remaining
   95%  (77/81)   17:09 elapsed   about 0:53 remaining
  100%  (81/81)   18:01 elapsed

  Best after zoom (.433, .203, inf), 2:12:01 total
    Lambda unit  0 --------------[===*===========]--- .8
    Lambda time  0 ------[==========*=====]---------- .4
    Lambda nn    0 [================================] .01

Visualisation of RMSE across lambda time and unit (X lowest RMSE, . highest RMSE)
     Each cell shows the lowest RMSE across all lambda_nn's searched in the adapted range.

Lambda time
     .278 | - - . - - - - - - |
          | + + + + + + + + + |
          | # # # # # # # # # |      +- Legend -----------------+
          | # # X # # # # # # |      |   | RMSE                 |
          | # # # # # # # # # |      | X | Best RMSE (7.08e-12) |
          | # # # # # # # # # |      | # | 7.08e-12 to 1.65e-11 |
          | # # # # # # # # # |      | + | 1.65e-11 to 3.54e-11 |
          | # # # # # # # # # |      | - | 3.54e-11 to 5.43e-11 |
    .0778 | # # # # # # # # # |      | . | 5.43e-11 to 6.38e-11 |
          +-------------------+      +--------------------------+
          .333            .733    Lambda unit

  Best over all evaluated lambdas (0, 0, .005)

Selected lambda: unit, time and nn = [0 ; 0 ; .005]
```

</details>

```
----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.04648
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |           1
-------------+--------------------------------------------------
 lambda_unit |      0.0000
 lambda_time |      0.0000
   lambda_nn |        .005
             |  (selected by loocv CV)
----------------------------------------------------------------
```
The adaptive search ends at the corner of the starting grid, (0, 0, .005), and the ATT moves further from zero, to 0.04648. The heat maps show why. The LOOCV RMSE is of the order 1e-12 to 1e-6 across the whole low-lambda region, so the criterion is flat there and cannot separate one triplet from another. With a single treated cell the adaptive refinement adds nothing over cycle.
### Block adoption
 
Under block adoption, units are treated simultaneously in the final 18 periods. Under `group(time)`, which is the default, this is one pooled spell of 270 cells (15 units × 18 periods):

```s
trop y unit time w_block, cv(resample, seed(1))
```

which returns 

```
Cross-validating lambdas using resample with cycle search, and 200 trials (seed 1).
To reduce resample computational time, reduce no of trials or set lambdas.
  cycle 1 of up to 50:  lambda_unit -> 1.2   lambda_time -> 1   lambda_nn -> .005   3:58 elapsed
  cycle 2 of up to 50:  lambda_unit -> .8   lambda_time -> 1   lambda_nn -> .1   17:03 elapsed
  cycle 3 of up to 50:  lambda_unit -> .8   lambda_time -> 1   lambda_nn -> .1   22:36 elapsed
  converged after 3 cycle(s), 22:36 total

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.00923
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          15
-------------+--------------------------------------------------
 lambda_unit |      0.8000
 lambda_time |      1.0000
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
  marginal:  lambda_time -> .2   lambda_unit -> 2   lambda_nn -> .005   0:34 elapsed
  cycle 1 of up to 50:  lambda_unit -> .2   lambda_time -> .2   lambda_nn -> .1   1:45 elapsed
  cycle 2 of up to 50:  lambda_unit -> .3   lambda_time -> .2   lambda_nn -> .25   2:50 elapsed
  cycle 3 of up to 50:  lambda_unit -> .3   lambda_time -> .2   lambda_nn -> .25   3:53 elapsed
  converged after 3 cycle(s), 3:53 total
Computing 270 group effects.
   35%  (95/270)   0:30 elapsed   about 0:55 remaining
   75%  (203/270)   1:02 elapsed   about 0:20 remaining
  100%  (270/270)   1:22 elapsed

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |    -0.00719
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          15
-------------+--------------------------------------------------
 lambda_unit |      0.3000
 lambda_time |      0.2000
   lambda_nn |         .25
             |  (selected by loocv CV)
----------------------------------------------------------------
(270 per-cell effects: e(group_grid) [unit x time], e(group_tau), e(group_info))
```
 
Both routes give an ATT close to zero, but the pooled estimate hides a lot of per-cell variation. We'll study heterogenous treatment effects more closely by using `e(group_grid)`:
 
```s
matrix list e(group_grid), format(%9.4f)
```
 
```
e(group_grid)[15,18]
          t31      t32      t33      t34      t35      t36      t37      t38      t39      t40      t41
 u97  -0.0365  -0.0446  -0.0600  -0.1164  -0.1181  -0.1522  -0.1930  -0.2088  -0.1897  -0.2381  -0.2843
 u98  -0.0157  -0.0556  -0.1148  -0.2970  -0.2135  -0.1952  -0.1934  -0.2359  -0.3000  -0.3242  -0.3739
 u99   0.1373   0.2239   0.3171   0.3889   0.4317   0.4656   0.4924   0.4637   0.3882   0.4053   0.4233
u100  -0.0467  -0.0649  -0.1164  -0.1465  -0.1191  -0.1114  -0.0971  -0.0906  -0.0496  -0.0071   0.0209
u101   0.0332   0.0430   0.0802   0.0764   0.0842   0.0618   0.0820   0.0906   0.1091   0.1328   0.1426
(output truncated)
```
 
The per-cell effects drift away from zero the further a cell sits from the onset of treatment, with opposite signs across units, so the pooled ATT of -0.007 is small only because they cancel. This is a mismatch between the LOOCV design and block adoption. LOOCV estimates placebo effects on control cells which are surrounded by untreated donor cells on both sides, so the criterion says nothing about extrapolating up to 18 periods past the last untreated period. Resample CV under `group(time)` suits block adoption, because the CV uses the actual treatment pattern on a subset of never-treated control units from the panel. As a result, we recommend using resample under block adoption with many treated periods. 
 
### Staggered adoption
 
Next, we study staggered adoption with three cohorts of five units. Under `group(time)`, the `detail` option prints a summary which includes which units are included in each cohort:
 
```s
trop y unit time w_stag, cv(resample, seed(1)) detail
```
 
```
Cross-validating lambdas using resample with cycle search, and 200 trials (seed 1).
To reduce resample computational time, reduce no of trials or set lambdas.
  cycle 1 of up to 50:  lambda_unit -> 1.2   lambda_time -> 1   lambda_nn -> .5   14:05 elapsed
  cycle 2 of up to 50:  lambda_unit -> 1.2   lambda_time -> 1   lambda_nn -> .5   24:42 elapsed
  converged after 2 cycle(s), 24:42 total

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.01913
             |  (no inference; vce(noinference))
-------------+--------------------------------------------------
     N units |         111
   T periods |          48
   N treated |          15
-------------+--------------------------------------------------
 lambda_unit |      1.2000
 lambda_time |      1.0000
   lambda_nn |          .5
             |  (selected by resample CV)
----------------------------------------------------------------
(3 per-spell effects in e(group_tau), e(group_weight), e(group_info), e(group_units))

time periods  start   end  units  cells        tau   cohort's units
      t21_48     21    48      5    140    -0.0185   97 98 99 100 101
      t31_48     31    48      5     90     0.1243   102 103 104 105 106
      t41_48     41    48      5     40    -0.0859   107 108 109 110 111
```
 
The per-cohort effects are also returned in `e(group_tau)` if you have not used the detail option (but the indexed cohort units are not):
 
```s
matrix list e(group_tau)
```
 
```
e(group_tau)[1,3]
        t21_48      t31_48      t41_48
r1  -.01845019   .12427219  -.08588071
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
         ATT |     0.00207
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
  100%  (200/200)   0:05 elapsed

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.02050
   Std. err. |     0.02094
      95% CI |   -0.02298    0.05610
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
  100%  (111/111)   0:02 elapsed

----------------------------------------------------------------
        TROP |  Triply Robust Panel estimator
-------------+--------------------------------------------------
         ATT |     0.02050
   Std. err. |     0.02245
      95% CI |   -0.02349    0.06450
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
