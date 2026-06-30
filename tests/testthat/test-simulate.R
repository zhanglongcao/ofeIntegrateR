test_that("simulate_ofe_trial returns expected structure and sizes", {
  sim <- simulate_ofe_trial(n_row = 10, n_col = 9, n_treat = 3,
                             treat_effects = c(0, 0.8, 1.6),
                             n_point_samples = 6, seed = 1)

  expect_named(sim, c("grid", "point_samples", "true_effects"))
  expect_equal(nrow(sim$grid), 10 * 9)
  expect_equal(nrow(sim$point_samples), 6)
  expect_true(all(c("row", "col", "treat", "point_true", "dense_response") %in%
                     names(sim$grid)))
  expect_true(all(c("row", "col", "point_true", "point_obs") %in%
                     names(sim$point_samples)))
  expect_equal(nlevels(sim$grid$treat), 3)
  expect_equal(unname(sim$true_effects), c(0, 0.8, 1.6))
})

test_that("simulate_ofe_trial errors on mismatched treat_effects length", {
  expect_error(
    simulate_ofe_trial(n_treat = 3, treat_effects = c(0, 1), seed = 1),
    "treat_effects"
  )
})

test_that("simulate_ofe_trial is reproducible with the same seed", {
  sim1 <- simulate_ofe_trial(n_row = 8, n_col = 6, n_point_samples = 5, seed = 42)
  sim2 <- simulate_ofe_trial(n_row = 8, n_col = 6, n_point_samples = 5, seed = 42)
  expect_equal(sim1$grid$dense_response, sim2$grid$dense_response)
  expect_equal(sim1$point_samples$point_obs, sim2$point_samples$point_obs)
})
