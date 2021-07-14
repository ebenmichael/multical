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

  # check if weights are post-stratification weights
  expect_equal(cal$weight, cal$target_count / cal$sample_count)

  imbal <- get_balance(cal, 2)$difference / sum(cal$target_count)
  expect_equal(imbal, numeric(length(imbal)))
})


test_that("Post stratification is exactfor 2nd order with cell level data", {

  cal <- multical(sample_count ~ X1 + X2, target_count, data_cell,
                  order = 2, lambda = 0, eps_rel = 1e-10, eps_abs = 1e-10)

  # check if weights are post-stratification weights
  expect_equal(cal$weight, cal$target_count / cal$sample_count)

  imbal <- get_balance(cal, 2)$difference  / sum(cal$target_count)
  expect_equal(imbal, numeric(length(imbal)))
})
