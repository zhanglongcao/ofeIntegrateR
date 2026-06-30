test_that("krige_point_samples returns predictions for all newdata rows", {
  sim <- simulate_ofe_trial(n_row = 12, n_col = 10, n_point_samples = 15, seed = 2)
  out <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")

  expect_equal(nrow(out), nrow(sim$grid))
  expect_true(all(c("point_obs_kriged", "point_obs_kriged_var") %in% names(out)))
  expect_false(anyNA(out$point_obs_kriged))
  expect_true(!is.null(attr(out, "variogram")))
})

test_that("krige_point_samples errors on missing columns", {
  sim <- simulate_ofe_trial(n_row = 10, n_col = 8, n_point_samples = 8, seed = 3)
  expect_error(
    krige_point_samples(sim$point_samples, sim$grid, value = "not_a_column"),
    "not found"
  )
  expect_error(
    krige_point_samples(sim$point_samples, sim$grid, value = "point_obs",
                         coords = c("x", "y")),
    "coords"
  )
})
