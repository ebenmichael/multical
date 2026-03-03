data(data_individual)
data(data_cell)

sample_ind <- data_individual %>% filter(response == 1)
pop_ind    <- data_individual  # all rows; target_count defaults to 1 per row

test_that("Calibration is exact with individual level data", {

  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                  order = 1, eps_rel = 1e-10, eps_abs = 1e-10)

  imbal <- get_balance(cal, 1)$difference /
    sum(cal$cells$target_count)
  expect_equal(imbal, numeric(length(imbal)))
})


test_that("Calibration is exact with cell level population data", {

  # sample is individual-level; population is the pre-aggregated cell table
  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, data_cell,
                  target_count = target_count,
                  order = 1, eps_rel = 1e-10, eps_abs = 1e-10)

  imbal <- get_balance(cal, 1)$difference /
    sum(cal$cells$target_count)
  expect_equal(imbal, numeric(length(imbal)))
})


test_that("Post stratification is exact for 2nd order with individual level data", {

  # With individual-level data every respondent has sample_count = 1, so the
  # regularization lambda * sum_i w_i^2 is symmetric within each (X1, X2) cell.
  # Any positive lambda uniquely determines the solution; by symmetry all
  # respondents in the same cell receive equal weights = PS weight.
  cal <- multical(~ X1 + X2, sample_ind, pop_ind,
                  order = 2, lambda = 0, eps_rel = 1e-10, eps_abs = 1e-10)

  cal_resp <- cal$cells %>%
    filter(sample_count != 0) %>%
    mutate(weight = weights(cal))

  ps_check <- cal_resp %>%
    group_by(X1, X2) %>%
    summarise(
      # within-cell weight SD should be ~0; coalesce for single-respondent cells
      weight_sd   = coalesce(sd(.data$weight), 0),
      mean_weight = mean(.data$weight),
      ps_weight   = sum(.data$target_count) / sum(.data$sample_count),
      .groups = "drop"
    )
  expect_equal(ps_check$weight_sd,   rep(0, nrow(ps_check)), tolerance = 1e-4)
  expect_equal(ps_check$mean_weight, ps_check$ps_weight,     tolerance = 1e-4)

  imbal <- get_balance(cal, 2)$difference /
    sum(cal$cells$target_count)
  expect_equal(imbal, numeric(length(imbal)), tolerance = 1e-4)
})


test_that("Post stratification is exact for 2nd order with cell level population data", {

  # Cell-level pop aggregated to X1 x X2; sample is individual-level respondents.
  # With sample_count = 1 for every respondent, the regularisation is symmetric
  # within cells, so each respondent in the same (X1, X2) cell gets the same
  # weight, which must equal pop_j / n_respondents_j.
  pop_2way <- data_cell %>%
    group_by(X1, X2) %>%
    summarise(target_count = sum(target_count), .groups = "drop")

  cal <- multical(~ X1 + X2, sample_ind, pop_2way,
                  target_count = target_count,
                  order = 2, lambda = 1e-8, eps_rel = 1e-10, eps_abs = 1e-10)

  n_resp_per_cell <- sample_ind %>%
    group_by(X1, X2) %>%
    summarise(n_resp = n(), .groups = "drop")

  cal_resp <- cal$cells %>%
    filter(sample_count != 0) %>%
    mutate(weight = weights(cal))

  ps_check <- cal_resp %>%
    group_by(X1, X2) %>%
    summarise(
      weight_sd   = coalesce(sd(.data$weight), 0),
      mean_weight = mean(.data$weight),
      .groups = "drop"
    ) %>%
    left_join(pop_2way,        by = c("X1", "X2")) %>%
    left_join(n_resp_per_cell, by = c("X1", "X2")) %>%
    mutate(ps_weight = target_count / n_resp)

  expect_equal(ps_check$weight_sd,   rep(0, nrow(ps_check)), tolerance = 1e-4)
  expect_equal(ps_check$mean_weight, ps_check$ps_weight,     tolerance = 1e-4)

  imbal <- get_balance(cal, 2)$difference /
    sum(cal$cells$target_count)
  expect_equal(imbal, numeric(length(imbal)), tolerance = 1e-4)
})


test_that("Individual-level and cell-level population data give the same weights", {

  # data_cell is the full 4-way aggregation of data_individual, so the two
  # population representations are equivalent and must yield identical weights.

  # order = 1 (raking): lambda is internally set to 0, fully determined
  cal_ind  <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                       order = 1, eps_rel = 1e-10, eps_abs = 1e-10)
  cal_cell <- multical(~ X1 + X2 + X3 + X4, sample_ind, data_cell,
                       target_count = target_count,
                       order = 1, eps_rel = 1e-10, eps_abs = 1e-10)

  expect_equal(weights(cal_ind), weights(cal_cell), tolerance = 1e-6)

  # order = 4, fixed lambda: higher-order post-stratification
  cal_ind4  <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                        order = 4, lambda = 1,
                        eps_rel = 1e-10, eps_abs = 1e-10)
  cal_cell4 <- multical(~ X1 + X2 + X3 + X4, sample_ind, data_cell,
                        target_count = target_count,
                        order = 4, lambda = 1,
                        eps_rel = 1e-10, eps_abs = 1e-10)

  expect_equal(weights(cal_ind4), weights(cal_cell4), tolerance = 1e-6)
})


test_that("Large lambda with base weights converges to base weights", {

  # With lambda -> Inf the penalty forces weights toward base_weights subject
  # to calibration constraints.  We first compute calibrated weights without
  # any base weights (these are constraint-feasible by construction), then use
  # those as base_weights.  At very large lambda the solution must be the
  # projection of base_weights onto the feasible set, which in this case is
  # base_weights itself, so the two runs should agree.
  cal0 <- multical(~ X1 + X2, sample_ind, pop_ind,
                   order = 2, lambda = 1e-8)

  base_wts <- weights(cal0)

  sample_ind_bw <- sample_ind %>% mutate(base_wt = base_wts)

  cal <- multical(~ X1 + X2, sample_ind_bw, pop_ind,
                  base_weights = base_wt,
                  order = 2, lambda = 1e6,
                  eps_rel = 1e-10, eps_abs = 1e-10)

  respondent_weights <- weights(cal)

  expect_equal(respondent_weights, base_wts, tolerance = 1e-2)
})


test_that("Base weights work correctly when pop-only cells exist", {

  # Drop all respondents from one level of X1 to create pop-only cells
  sample_partial <- sample_ind %>% filter(!(X1 == 1 & X2 == 1))
  base_wts <- rep(1, nrow(sample_partial))
  sample_partial_bw <- sample_partial %>% mutate(base_wt = base_wts)

  # Should run without error and produce finite weights for all respondents
  expect_no_error(
    cal <- multical(~ X1 + X2, sample_partial_bw, pop_ind,
                    base_weights = base_wt,
                    order = 2, lambda = 0,
                    eps_rel = 1e-10, eps_abs = 1e-10)
  )

  respondent_weights <- weights(cal)

  expect_true(all(is.finite(respondent_weights)))
})


test_that("multical runs without error when lambda is NULL", {

  # Default behaviour: fit over a grid of lambda values
  expect_no_error(
    cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                    order = 4, eps_rel = 1e-8, eps_abs = 1e-8)
  )

  # Should return a multical object with multiple lambdas and a weight matrix
  expect_s3_class(cal, "multical")
  expect_gt(length(cal$lambda), 1L)
  expect_equal(nrow(cal$weights), nrow(sample_ind))
  expect_equal(ncol(cal$weights), length(cal$lambda))

  # weights() should return a vector of length equal to n_respondents
  w <- weights(cal)
  expect_length(w, nrow(sample_ind))
  expect_true(all(is.finite(w)))
})


test_that("weights() returns weights at the correct default lambda", {

  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                  order = 4, eps_rel = 1e-8, eps_abs = 1e-8)

  idx <- cal$default_lambda_idx

  # weights() must exactly match the column of the weight matrix at default_lambda_idx
  expect_equal(weights(cal), cal$weights[, idx])

  # default_lambda_idx must reproduce what select_default_lambda returns when
  # called directly with the same inputs
  expect_equal(
    idx,
    multical:::select_default_lambda(
      cal$weights, cal$cells, cal$lambda, cal$order, cal$balance_threshold
    )
  )

  # The default lambda should not be the worst-balance one (index 1, largest lambda),
  # since it should achieve at least balance_threshold of the possible gain
  expect_gt(idx, 1L)

  # Weights should be positive and finite
  w <- weights(cal)
  expect_true(all(w >= 0))
  expect_true(all(is.finite(w)))
})


test_that("estimate() returns correct point estimate and SE for Hajek method", {

  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                  order = 1, eps_rel = 1e-10, eps_abs = 1e-10)

  # --- via bare column name + data argument ---
  result_data <- estimate(cal, y, data = sample_ind, method = "hajek")

  w   <- weights(cal)
  est <- sum(w * sample_ind$y) / sum(w)
  se  <- sqrt(sum(w ^ 2 * (sample_ind$y - est) ^ 2)) / sum(w)

  expect_equal(result_data$estimate, est, tolerance = 1e-10)
  expect_equal(result_data$se, se, tolerance = 1e-10)
  expect_equal(result_data$method, "hajek")
  expect_equal(result_data$lambda, cal$lambda[cal$default_lambda_idx])
  expect_equal(nrow(result_data), 1L)

  # --- via inline vector expression (no data argument) ---
  result_vec <- estimate(cal, sample_ind$y, method = "hajek")
  expect_equal(result_vec$estimate, est, tolerance = 1e-10)

  # SE should be positive and finite
  expect_gt(result_data$se, 0)
  expect_true(is.finite(result_data$se))
})


test_that("estimate() respects lambda_idx and balance_threshold overrides", {

  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                  order = 4, eps_rel = 1e-8, eps_abs = 1e-8)

  # explicit lambda_idx
  idx <- 5L
  result_idx <- estimate(cal, sample_ind$y, lambda_idx = idx)
  w_idx <- weights(cal, lambda_idx = idx)
  expect_equal(result_idx$estimate,
               sum(w_idx * sample_ind$y) / sum(w_idx), tolerance = 1e-10)
  expect_equal(result_idx$lambda, cal$lambda[idx])

  # balance_threshold override
  result_bt <- estimate(cal, sample_ind$y, balance_threshold = 0.5)
  idx_bt <- select_default_lambda(
    cal$weights, cal$cells, cal$lambda, cal$order, 0.5
  )
  w_bt <- weights(cal, lambda_idx = idx_bt)
  expect_equal(result_bt$estimate,
               sum(w_bt * sample_ind$y) / sum(w_bt), tolerance = 1e-10)
})


test_that("estimate() linearized method gives correct point estimate and smaller SE when y is correlated with covariates", {

  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                  order = 4)

  result_linearized <- estimate(cal, y, data = sample_ind, method = "linearized",
                                order = 4)

  # Point estimate is still the weighted mean
  w   <- weights(cal)
  est <- sum(w * sample_ind$y) / sum(w)
  expect_equal(result_linearized$estimate, est, tolerance = 1e-10)
  expect_equal(result_linearized$method, "linearized")
  expect_equal(result_linearized$lambda, cal$lambda[cal$default_lambda_idx])

  # SE must be positive and finite
  expect_gt(result_linearized$se, 0)
  expect_true(is.finite(result_linearized$se))

  # Manually verify SE using regression residuals at reg_order = order = 4
  cov_cols   <- setdiff(colnames(cal$cells),
                        c("sample_count", "target_count", "base_weight"))
  cells_resp <- cal$cells[cal$cells$sample_count != 0, cov_cols, drop = FALSE]
  X      <- model.matrix(~ .^4, data = cells_resp)
  resids <- lm.fit(X, sample_ind$y)$residuals
  se_manual <- sqrt(sum(w ^ 2 * resids ^ 2)) / sum(w)
  expect_equal(result_linearized$se, se_manual, tolerance = 1e-10)

  # Linearized SE should be <= Hajek SE when y is correlated with covariates
  result_hajek <- estimate(cal, y, data = sample_ind, method = "hajek")
  expect_lte(result_linearized$se, result_hajek$se)
})


test_that("estimate() linearized method with use_ridge uses glmnet CV residuals", {

  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                  order = 1, eps_rel = 1e-10, eps_abs = 1e-10)
  w <- weights(cal)

  # use_ridge = FALSE (default) should match plain OLS
  result_ols   <- estimate(cal, y, data = sample_ind, method = "linearized")
  result_noridge <- estimate(cal, y, data = sample_ind,
                             method = "linearized", use_ridge = FALSE)
  expect_equal(result_noridge$se, result_ols$se, tolerance = 1e-10)

  # use_ridge = TRUE: verify manually using cv.glmnet
  cov_cols   <- setdiff(colnames(cal$cells),
                        c("sample_count", "target_count", "base_weight"))
  cells_resp <- cal$cells[cal$cells$sample_count != 0, cov_cols, drop = FALSE]
  X        <- model.matrix(~ ., data = cells_resp)
  X_noInt  <- X[, colnames(X) != "(Intercept)", drop = FALSE]

  set.seed(1011)
  cv_fit   <- glmnet::cv.glmnet(X_noInt, sample_ind$y, alpha = 0)
  fitted   <- as.numeric(predict(cv_fit, newx = X_noInt, s = "lambda.min"))
  resids   <- sample_ind$y - fitted
  se_manual <- sqrt(sum(w ^ 2 * resids ^ 2)) / sum(w)

  set.seed(1011)
  result_ridge <- estimate(cal, y, data = sample_ind,
                           method = "linearized", use_ridge = TRUE)

  expect_equal(result_ridge$se, se_manual, tolerance = 1e-10)
  expect_equal(result_ridge$estimate,
               sum(w * sample_ind$y) / sum(w), tolerance = 1e-10)
  expect_gt(result_ridge$se, 0)
  expect_true(is.finite(result_ridge$se))
})


test_that("estimate() greg with order = 1 gives the same point estimate as linearized (with and without ridge)", {
  # When order-1 calibration constraints are satisfied, the GREG correction
  # term (X_target_bar - X_sample_bar) %*% beta is zero, so the greg point
  # estimate collapses to the Hajek weighted mean -- same as linearized.
  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                  order = 1, eps_rel = 1e-10, eps_abs = 1e-10)

  result_lin  <- estimate(cal, y, data = sample_ind, method = "linearized")
  result_greg <- estimate(cal, y, data = sample_ind, method = "greg", order = 1)

  expect_equal(result_greg$estimate, result_lin$estimate, tolerance = 1e-6)

  # Same should hold when ridge is used
  set.seed(1011)
  result_greg_ridge <- estimate(cal, y, data = sample_ind,
                                method = "greg", order = 1, use_ridge = TRUE)
  expect_equal(result_greg_ridge$estimate, result_lin$estimate, tolerance = 1e-6)
})


test_that("estimate() drp method returns correctly structured output", {
  skip_if_not_installed("xgboost")
  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind, order = 1)

  set.seed(42)
  result_drp <- estimate(cal, y, data = sample_ind, method = "drp",
                         nrounds = 50)

  expect_s3_class(result_drp, "data.frame")
  expect_equal(nrow(result_drp), 1)
  expect_named(result_drp, c("estimate", "se", "lambda", "method"))
  expect_true(is.finite(result_drp$estimate))
  expect_gt(result_drp$se, 0)
  expect_true(is.finite(result_drp$se))
  expect_equal(result_drp$method, "drp")
  expect_equal(result_drp$lambda, cal$lambda[cal$default_lambda_idx])
})


test_that("estimate() drp method is reproducible with set.seed()", {
  skip_if_not_installed("xgboost")
  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind, order = 1)

  set.seed(123)
  result1 <- estimate(cal, y, data = sample_ind, method = "drp", nrounds = 50)
  set.seed(123)
  result2 <- estimate(cal, y, data = sample_ind, method = "drp", nrounds = 50)

  expect_equal(result1$estimate, result2$estimate)
  expect_equal(result1$se, result2$se)
})


test_that("estimate() handles factor y by estimating each level as a binary outcome", {
  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind, order = 1)

  sample_ind$y_fac <- factor(sample_ind$y, levels = c(0, 1))
  result <- estimate(cal, y_fac, data = sample_ind)

  y_fac <- sample_ind$y_fac
  lvls <- levels(y_fac)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), length(lvls))
  expect_named(result, c("level", "estimate", "se", "lambda", "method"))
  expect_equal(result$level, lvls)

  # each row should match a direct binary-indicator call
  for (i in seq_along(lvls)) {
    y_bin    <- as.numeric(y_fac == lvls[i])
    expected <- estimate(cal, y_bin)
    expect_equal(result$estimate[i], expected$estimate, tolerance = 1e-10)
    expect_equal(result$se[i],       expected$se,       tolerance = 1e-10)
  }
})


test_that("estimate() coerces character y to factor and gives the same result", {
  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind, order = 1)
  sample_ind$y_fac <- factor(sample_ind$y, levels = c(0, 1))
  y_chr <- as.character(sample_ind$y_fac)

  result_fac <- estimate(cal, y_fac, data = sample_ind)
  result_chr <- estimate(cal, y_chr, data = sample_ind)

  expect_equal(result_fac, result_chr)
})


test_that("estimate() errors for non-numeric, non-factor, non-character y", {
  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind, order = 1)
  expect_error(estimate(cal, as.list(sample_ind$y)), "`y` must evaluate to")
})
