---
title: "multical"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{multical}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




# `multical`: Multilevel calibration weighting for surveys

## Installation
You can install `multical` from github using `devtools`.


```r
## Install devtools if not already installed
install.packages("devtools", repos='http://cran.us.r-project.org')
## Install multical from github
devtools::install_github("ebenmichael/multical")
```


## Data

`multical` takes two separate data frames as inputs:

- **`sample_data`**: individual-level data for survey respondents. One row per respondent.
- **`pop_data`**: data describing the target population. This can be individual-level (one row per population member) or pre-aggregated to the cell level (one row per unique combination of covariate values, with a column giving the count per cell).

We'll use two contrived examples shipped with the package: `data_individual` has individual-level information on covariates, response, and outcome; `data_cell` has cell-level information.


```r
library(multical)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
data(data_individual)
data(data_cell)
```

`data_individual` records 4 covariates, whether each individual responded (`response`) to the survey, whether they are in the target population (`intarget`), and their answer to a binary question `y`.


```r
head(data_individual)
#>   X1 X2 X3 X4 response intarget y
#> 1  3  2  2  4        1        1 0
#> 2  2  2  2  4        1        1 0
#> 3  1  1  2  1        1        1 1
#> 4  3  1  1  4        1        1 0
#> 5  2  2  1  1        1        1 0
#> 6  1  1  1  1        1        1 1
```

`data_cell` records each unique combination of the 4 covariates, the number of respondents in that cell (`sample_count`), the number of individuals in the target population in that cell (`target_count`), and the mean outcome in the cell.


```r
head(data_cell)
#> # A tibble: 6 × 7
#>   X1    X2    X3    X4    sample_count target_count     y
#>   <fct> <fct> <fct> <fct>        <dbl>        <dbl> <dbl>
#> 1 1     1     1     1                8           50 0.42 
#> 2 1     1     1     2                5           21 0.810
#> 3 1     1     1     3                2           17 0.765
#> 4 1     1     1     4               33          281 0.530
#> 5 1     1     2     1               10           29 0.552
#> 6 1     1     2     2                9           23 0.565
```

For the examples below we'll split `data_individual` into respondents (our sample) and the full population.


```r
sample_ind <- data_individual %>% filter(response == 1)
pop_ind    <- data_individual   # all rows; each counts as 1 population member
```


## Getting multilevel calibration weights

The main function is `multical`. Its first three arguments are:

1. `formula`: a **right-hand-side only** formula (`~ covariates`) specifying which variables to calibrate on
2. `sample_data`: the respondent data frame
3. `pop_data`: the population data frame (individual or cell level)

An optional `target_count` argument names a column in `pop_data` giving the count or weight per row. When omitted, each row of `pop_data` is treated as one population member.

### Raking on first-order margins (individual-level population)

Setting `order = 1` calibrates only on the main-effect margins of each covariate.


```r
out <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind, order = 1)
#> Error in terms(Formula::Formula(formula), rhs = 1)[[3]]: subscript out of bounds
out
#> Error in eval(expr, envir, enclos): object 'out' not found
```

`multical` returns a `multical` object. The key fields are:

- `out$weights`: matrix of calibration weights (one column per lambda value, one row per respondent, in the same order as `sample_data`)
- `out$cells`: data frame of covariate values plus `sample_count`, `target_count`, and `base_weight` for every unit (respondents and any pop-only cells)
- `out$lambda`: the lambda values used

Call `weights(out)` to extract respondent weights at the auto-selected default lambda.

We can inspect marginal balance using `get_balance`:


```r
get_balance(out, 1)
#> Error in eval(expr, envir, enclos): object 'out' not found
```

To estimate the population mean of an outcome, use the `estimate()` function
(defaults to the linearized estimator — see the [Estimating population means](#estimating-population-means-with-estimate) section below for all options):


```r
estimate(out, y, data = sample_ind)
#> Error in estimate(out, y, data = sample_ind): could not find function "estimate"
```

### Higher-order calibration with a fixed lambda

Setting `order > 1` includes interaction terms up to that order. The parameter `lambda` controls the trade-off between exact post-stratification (small `lambda`) and lower variance (large `lambda`).


```r
out <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind, order = 4, lambda = 1)
#> Error in terms(Formula::Formula(formula), rhs = 1)[[3]]: subscript out of bounds

rbind(head(get_balance(out, 4)), tail(get_balance(out, 4)))
#> Error in eval(expr, envir, enclos): object 'out' not found
```

### Selecting lambda automatically

By default (no `lambda` supplied), `multical` solves over a grid of values from `lambda_max` down to `lambda_max * lambda_min_ratio`. Use `get_balance_v_sample_size` to trace the balance-effective-sample-size trade-off and choose an appropriate value.


```r
out <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind)
#> Error in terms(Formula::Formula(formula), rhs = 1)[[3]]: subscript out of bounds

# plot() shows the balance vs. effective-sample-size trade-off and
# highlights the auto-selected default lambda in red
plot(out)
#> Error in eval(expr, envir, enclos): object 'out' not found
```

By default, the `weights()` function will extract the weights corresponding to the value of `lambda` that achieves 90% of the total possible balance gain. You can override this by either (i) supplying a specific index to choose from (e.g. `weights(out, lambda_idx = 2)` returns weights corresponding to `out$lambda[2]`) or (ii) supplying a specific amount of balance gain to target (e.g. `weights(out, balance_gain_target = 0.95)` returns weights that give at least 95% of the total possible balance gain).

For large datasets with many covariates, start with a low `order` and increase as needed — the number of interaction terms grows quickly.

### Using pre-aggregated cell-level population data

If the population is available as a pre-aggregated cell table (e.g. from a census), pass it as `pop_data` and name the count column via `target_count`. The sample is still individual-level.


```r
out <- multical(~ X1 + X2 + X3 + X4, sample_ind, data_cell,
                target_count = target_count,
                order = 4, lambda = 1)
#> Error in terms(Formula::Formula(formula), rhs = 1)[[3]]: subscript out of bounds
```

Estimation works the same way (using the default linearized estimator):


```r
estimate(out, y, data = sample_ind)
#> Error in estimate(out, y, data = sample_ind): could not find function "estimate"
```


## Estimating population means with `estimate()`

The `estimate()` function computes a population mean estimate and a
linearization-based standard error from a fitted `multical` object. The
`method` argument selects the estimator:

| Method | Description |
|--------|-------------|
| `"linearized"` (default) | Point estimate is the weighted mean. Standard errors are the residuals of an OLS or ridge regression of the outcome on interactions up to the specified order. |
| `"hajek"` | Point estimate is the weighted mean. Standard error is the standard sandwich standard error.|
| `"greg"` | GREG estimator; uses OLS/ridge regression for model assistance to reduce bias. Standard errors are based off the residuals of the regression |
| `"drp"` | Double regression and post-stratification; the same implementation as `"greg"` but uses cross-fitted gradient-boosted trees (xgboost) for the outcome model. |

All methods return a one-row data frame with columns `estimate`, `se`, `lambda`,
and `method`.


```r
out <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind, order = 4)
#> Error in terms(Formula::Formula(formula), rhs = 1)[[3]]: subscript out of bounds
```

### Linearized estimator (default)

The default. Fits an OLS model of `y` on the covariate design matrix and uses
the residuals in the SE formula.
The `order` argument sets the interaction order in the regression (defaults to 1).


```r
estimate(out, y, data = sample_ind)
#> Error in estimate(out, y, data = sample_ind): could not find function "estimate"
```

You can set the interaction order explicitly:


```r
estimate(out, y, data = sample_ind, method = "linearized", order = 2)
#> Error in estimate(out, y, data = sample_ind, method = "linearized", order = 2): could not find function "estimate"
```

Pass `use_ridge = TRUE` to fit the outcome model with ridge regression instead
(penalty chosen by 10-fold cross-validation via `glmnet`):


```r
estimate(out, y, data = sample_ind, method = "linearized", order = 2,
         use_ridge = TRUE)
#> Error in estimate(out, y, data = sample_ind, method = "linearized", order = 2, : could not find function "estimate"
```

### Hajek estimator

Returns the calibration-weighted mean of `y` and a  standard sandwich SE,
without any regression adjustment.


```r
estimate(out, y, data = sample_ind, method = "hajek")
#> Error in estimate(out, y, data = sample_ind, method = "hajek"): could not find function "estimate"
```

### GREG estimator

The GREG estimator uses regression predictions on the target population to form
the point estimate. The `order` and `use_ridge` arguments work the same way as for
`"linearized"`.


```r
estimate(out, y, data = sample_ind, method = "greg", order = 2)
#> Error in estimate(out, y, data = sample_ind, method = "greg", order = 2): could not find function "estimate"
```

### DRP estimator

The double regression and post-stratification (DRP) estimator replaces the parametric regression
with cross-fitted gradient-boosted trees (xgboost). Additional arguments
(e.g. `nrounds`, `eta`) are forwarded directly to `xgboost`.


```r
estimate(out, y, data = sample_ind, method = "drp", nrounds = 200)
#> Error in estimate(out, y, data = sample_ind, method = "drp", nrounds = 200): could not find function "estimate"
```

### Selecting lambda

By default `estimate()` uses the auto-selected default lambda (the one that
achieves at least 90% of the possible balance gain). You can override this with
`lambda_idx` (an integer index into `out$lambda`) or `balance_threshold` (a
value in (0, 1) that re-runs lambda selection):


```r
# use the weights at index 5 in the lambda grid
estimate(out, y, data = sample_ind, lambda_idx = 5)
#> Error in estimate(out, y, data = sample_ind, lambda_idx = 5): could not find function "estimate"

# re-select lambda targeting 95% balance gain
estimate(out, y, data = sample_ind, balance_threshold = 0.95)
#> Error in estimate(out, y, data = sample_ind, balance_threshold = 0.95): could not find function "estimate"
```
