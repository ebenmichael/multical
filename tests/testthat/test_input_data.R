data(data_individual)

test_that("multical runs with grouped sample data", {

  sample_grp <- data_individual %>%
    filter(response == 1) %>%
    group_by(across(contains("X")))

  pop_ind <- data_individual

  expect_error(multical(~ X1 + X2 + X3 + X4, sample_grp, pop_ind,
                        order = 1, eps_rel = 1e-10, eps_abs = 1e-10),
               NA)
})

test_that("multical runs with grouped population data", {

  sample_ind <- data_individual %>% filter(response == 1)

  pop_grp <- data_individual %>% group_by(across(contains("X")))

  expect_error(multical(~ X1 + X2 + X3 + X4, sample_ind, pop_grp,
                        order = 1, eps_rel = 1e-10, eps_abs = 1e-10),
               NA)
})

# ---------------------------------------------------------------------------
# estimate.multical subset parameter tests
# ---------------------------------------------------------------------------

# Shared calibration object used across several estimate subset tests
local({
  sample_ind <- data_individual %>% filter(response == 1)
  pop_ind    <- data_individual
  cal        <<- multical(~ X1 + X2 + X3 + X4, sample_ind, pop_ind,
                          order = 1, eps_rel = 1e-10, eps_abs = 1e-10)
  n_resp     <<- cal$n_respondents
  set.seed(42)
  y_full     <<- rnorm(n_resp)
})

test_that("subset = NULL gives same result as no subset argument", {
  res_default <- estimate(cal, y_full, method = "hajek")
  res_null    <- estimate(cal, y_full, method = "hajek", subset = NULL)
  expect_equal(res_default, res_null)
})

test_that("subset via logical mask returns correct estimate", {
  mask   <- seq_len(n_resp) %% 2 == 1   # odd-indexed respondents
  y_sub  <- y_full[mask]
  w_sub  <- weights(cal)[mask]
  expected_est <- sum(w_sub * y_sub) / sum(w_sub)

  res <- estimate(cal, y_full, method = "hajek", subset = mask)
  expect_equal(res$estimate, expected_est, tolerance = 1e-10)
})

test_that("subset via integer index vector returns correct estimate", {
  idx    <- c(1L, 3L, 5L, 7L, 9L)
  y_sub  <- y_full[idx]
  w_sub  <- weights(cal)[idx]
  expected_est <- sum(w_sub * y_sub) / sum(w_sub)

  res <- estimate(cal, y_full, method = "hajek", subset = idx)
  expect_equal(res$estimate, expected_est, tolerance = 1e-10)
})

test_that("NA in y outside subset does not warn and result is unchanged", {
  y_na       <- y_full
  mask       <- seq_len(n_resp) %% 2 == 1
  # Place NA at an even position (NOT in the subset)
  na_pos     <- which(!mask)[1]
  y_na[na_pos] <- NA

  expect_no_warning(
    res <- estimate(cal, y_na, method = "hajek", subset = mask)
  )
  res_clean <- estimate(cal, y_full, method = "hajek", subset = mask)
  expect_equal(res$estimate, res_clean$estimate, tolerance = 1e-10)
})

test_that("NA in y inside subset warns and drops those rows", {
  y_na   <- y_full
  mask   <- seq_len(n_resp) %% 2 == 1
  na_pos <- which(mask)[1]   # first retained position
  y_na[na_pos] <- NA

  expect_warning(
    res <- estimate(cal, y_na, method = "hajek", subset = mask),
    "1 observation\\(s\\) with missing outcome \\(NA\\) excluded from estimation\\."
  )

  # Result should match estimating on the subset with na_pos also excluded
  keep_manual        <- mask
  keep_manual[na_pos] <- FALSE
  y_sub  <- y_full[keep_manual]
  w_sub  <- weights(cal)[keep_manual]
  expected_est <- sum(w_sub * y_sub) / sum(w_sub)
  expect_equal(res$estimate, expected_est, tolerance = 1e-10)
})

test_that("factor y with NA: warning fires once, not once per level", {
  y_fac  <- factor(sample(c("a", "b", "c"), n_resp, replace = TRUE))
  mask   <- seq_len(n_resp) %% 2 == 1
  na_pos <- which(mask)[1]
  y_fac[na_pos] <- NA

  # Expect exactly one warning (not three, one per level)
  warns <- character(0)
  withCallingHandlers(
    estimate(cal, y_fac, method = "hajek", subset = mask),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(length(warns), 1L)
  expect_match(warns[1], "excluded from estimation")
})

test_that("invalid subset type stops with error", {
  expect_error(estimate(cal, y_full, method = "hajek", subset = "bad"),
               "`subset` must be NULL")
})

test_that("logical subset of wrong length stops with error", {
  expect_error(estimate(cal, y_full, method = "hajek",
                        subset = rep(TRUE, n_resp + 1)),
               "`subset` must have length equal to the number of respondents")
})

test_that("integer subset with out-of-range index stops with error", {
  expect_error(estimate(cal, y_full, method = "hajek",
                        subset = c(1L, n_resp + 1L)),
               "`subset` contains out-of-range indices")
})