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