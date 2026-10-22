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

test_that("krige_point_samples falls back when gstat returns a negative range", {
  # gstat::fit.variogram() can return a "singular model" with a negative
  # range on clustered samples; gstat::krige() then stops. The constrained
  # fit must be used instead.
  sim <- simulate_ofe_trial(n_row = 12, n_col = 10, n_point_samples = 15, seed = 2)
  bad_fit <- function(object, model, ...) {
    warning("singular model in variogram fit")
    m <- gstat::vgm(psill = 1, model = "Exp", range = 1, nugget = 0.5)
    m$range[m$model == "Exp"] <- -0.2   # what gstat returned on a real sample
    m
  }
  local_mocked_bindings(fit.variogram = bad_fit, .package = "gstat")
  out <- expect_no_error(
    krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
  )
  vg <- attr(out, "variogram")
  expect_true(all(vg$range[vg$model != "Nug"] > 0))
  expect_equal(attr(out, "variogram_source"), "ofeIntegrateR constrained fit")
  expect_false(anyNA(out$point_obs_kriged))
})
