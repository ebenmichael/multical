data(data_individual)
data(data_cell)

sample_ind <- data_individual %>% filter(response == 1)
pop_ind    <- data_individual  # all rows; target_count defaults to 1 per row

test_that("Calibration is exact with individual level data", {

  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                  order = 1, eps_rel = 1e-10, eps_abs = 1e-10)

  imbal <- get_balance(cal, 1)$difference /
    sum(cal$target_count[cal$lambda == cal$lambda[1]])
  expect_equal(imbal, numeric(length(imbal)))
})


test_that("Calibration is exact with cell level population data", {

  # sample is individual-level; population is the pre-aggregated cell table
  cal <- multical(~ X1 + X2 + X3 + X4, sample_ind, data_cell,
                  target_count = target_count,
                  order = 1, eps_rel = 1e-10, eps_abs = 1e-10)

  imbal <- get_balance(cal, 1)$difference /
    sum(cal$target_count[cal$lambda == cal$lambda[1]])
  expect_equal(imbal, numeric(length(imbal)))
})


test_that("Post stratification is exactfor 2nd order with individual level data", {

  # With individual-level data every respondent has sample_count = 1, so the
  # regularization lambda * sum_i w_i^2 is symmetric within each (X1, X2) cell.
  # Any positive lambda uniquely determines the solution; by symmetry all
  # respondents in the same cell receive equal weights = PS weight.
  cal <- multical(~ X1 + X2, sample_ind, pop_ind,
                  order = 2, lambda = 0, eps_rel = 1e-10, eps_abs = 1e-10)

  ps_check <- cal %>%
    filter(.data$sample_count != 0) %>%
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
    sum(cal$target_count[cal$lambda == cal$lambda[1]])
  expect_equal(imbal, numeric(length(imbal)), tolerance = 1e-4)
})


test_that("Post stratification is exactfor 2nd order with cell level population data", {

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

  ps_check <- cal %>%
    filter(.data$sample_count != 0) %>%
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
    sum(cal$target_count[cal$lambda == cal$lambda[1]])
  expect_equal(imbal, numeric(length(imbal)), tolerance = 1e-4)
})
