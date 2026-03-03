context("Test that balancing is working")

data(data_individual)
data(data_cell)

test_that("Calibration is exact with individual level data", {

  cal <- multical(response ~ X1 + X2 + X3 + X4, 1 - response, data_individual,
                  order = 1, eps_rel = 1e-10, eps_abs = 1e-10)
  
  imbal <- get_balance(cal, 1)$difference / sum(cal$target_count)
  expect_equal(imbal, numeric(length(imbal)))
})


test_that("Calibration is exact with cell level data", {

  cal <- multical(sample_count ~ X1 + X2 + X3 + X4, target_count, data_cell,
                  order = 1, eps_rel = 1e-10, eps_abs = 1e-10)
  
  imbal <- get_balance(cal, 1)$difference / sum(cal$target_count)
  expect_equal(imbal, numeric(length(imbal)))
})


test_that("Post stratification is exactfor 2nd order with individual level data", {

  cal <- multical(response ~ X1 + X2, 1 - response, data_individual,
                  order = 2, lambda = 0, eps_rel = 1e-10, eps_abs = 1e-10)

  ps_check <- cal %>%
    group_by(X1, X2) %>%
    summarise(
      # within-cell weight SD should be ~0 (all respondents get the same weight);
      # coalesce to 0 for any single-respondent cell where sd() returns NA
      weight_sd  = coalesce(sd(.data$weight[.data$sample_count != 0]), 0),
      # mean respondent weight should equal the PS ratio computed from all rows
      mean_weight = mean(.data$weight[.data$sample_count != 0]),
      ps_weight   = sum(.data$target_count) / sum(.data$sample_count),
      .groups = "drop"
    )
  expect_equal(ps_check$weight_sd,   rep(0, nrow(ps_check)), tolerance = 1e-4)
  expect_equal(ps_check$mean_weight, ps_check$ps_weight,     tolerance = 1e-4)

  imbal <- get_balance(cal, 2)$difference /
    sum(cal$target_count[cal$lambda == cal$lambda[1]])
  expect_equal(imbal, numeric(length(imbal)), tolerance = 1e-4)
})


test_that("Post stratification is exactfor 2nd order with cell level data", {

  # data_cell is at 4-covariate granularity, so multiple rows share the same
  # (X1, X2) but have *different* sample_count values.  The regularisation
  # penalty lambda * sample_count_i * w_i^2 is therefore not symmetric within
  # groups, so individual-row weights are not forced equal and don't match the
  # PS ratio.  The correct check is at the aggregate level.
  cal <- multical(sample_count ~ X1 + X2, target_count, data_cell,
                  order = 2, lambda = 1e-8, eps_rel = 1e-10, eps_abs = 1e-10)

  ps_check <- cal %>%
    group_by(X1, X2) %>%
    summarise(weighted_sample = sum(.data$weight * .data$sample_count),
              target_total    = sum(.data$target_count),
              .groups = "drop")
  expect_equal(ps_check$weighted_sample, ps_check$target_total, tolerance = 1e-4)

  imbal <- get_balance(cal, 2)$difference /
    sum(cal$target_count[cal$lambda == cal$lambda[1]])
  expect_equal(imbal, numeric(length(imbal)), tolerance = 1e-4)
})
