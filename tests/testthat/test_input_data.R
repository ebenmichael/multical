context("Test that multical handles input data well")

data(data_individual)

test_that("multical runs with grouped data", {

  tmpdat <- data_individual %>% group_by(across(contains("X"))) %>%
    summarise(sample_count = sum(response), target_count = sum(intarget),
              y = mean(y))

  expect_error(multical(sample_count ~ X1 + X2 + X3 + X4, target_count, tmpdat, 
                        order = 1, eps_rel = 1e-10, eps_abs = 1e-10),
               NA)

})