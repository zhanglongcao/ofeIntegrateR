test_that("fit_integrated_kriged (engine = 'lm') returns sensible treatment contrasts", {
  sim <- simulate_ofe_trial(n_row = 30, n_col = 21, n_treat = 3,
                             treat_effects = c(0, 0.8, 1.6),
                             n_point_samples = 25, seed = 11)
  krieged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")

  fit <- fit_integrated_kriged(krieged, response = "dense_response",
                                treat = "treat", covariate = "point_obs_kriged",
                                row = "row", col = "col", engine = "lm")
  expect_s3_class(fit, "lm")

  fx <- extract_fixed_effects(fit)
  expect_true(all(c("term", "estimate", "se") %in% names(fx)))

  # True effects are A = 0 < B = 0.8 < C = 1.6. With contiguous treatment
  # strips confounded against a smooth spatial covariate, exact recovery
  # from a single simulated draw is not guaranteed, but the estimates
  # should be finite, positive, and ordered B < C.
  b_est <- fx$estimate[fx$term == "treatB"]
  c_est <- fx$estimate[fx$term == "treatC"]
  expect_true(is.finite(b_est) && is.finite(c_est))
  expect_gt(b_est, 0)
  expect_gt(c_est, b_est)
})

test_that("fit_integrated_kriged errors on missing columns", {
  sim <- simulate_ofe_trial(n_row = 10, n_col = 8, n_point_samples = 6, seed = 4)
  expect_error(
    fit_integrated_kriged(sim$grid, response = "not_a_column",
                           treat = "treat", covariate = "point_true",
                           engine = "lm"),
    "not found"
  )
})

test_that("fit_integrated_kriged (engine = 'asreml') requires the asreml package", {
  skip_if(requireNamespace("asreml", quietly = TRUE),
          "asreml is installed; informative-error path not applicable")
  sim <- simulate_ofe_trial(n_row = 10, n_col = 8, n_point_samples = 6, seed = 5)
  krieged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
  expect_error(
    fit_integrated_kriged(krieged, response = "dense_response",
                           treat = "treat", covariate = "point_obs_kriged",
                           engine = "asreml"),
    "asreml"
  )
})

test_that("fit_integrated_kriged (engine = 'asreml') runs when asreml is available", {
  skip_if_not_installed("asreml")
  sim <- simulate_ofe_trial(n_row = 30, n_col = 21, n_treat = 3,
                             treat_effects = c(0, 0.8, 1.6),
                             n_point_samples = 25, seed = 11)
  krieged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
  fit <- fit_integrated_kriged(krieged, response = "dense_response",
                                treat = "treat", covariate = "point_obs_kriged",
                                row = "row", col = "col", engine = "asreml")
  expect_s3_class(fit, "asreml")
  fx <- extract_fixed_effects(fit)
  expect_true(all(c("term", "estimate", "se") %in% names(fx)))

  # asreml labels factor-level coefficients "treat_B" (not lm's "treatB")
  b_est <- fx$estimate[grepl("treat_B$", fx$term)]
  c_est <- fx$estimate[grepl("treat_C$", fx$term)]
  expect_length(b_est, 1)
  expect_length(c_est, 1)
  expect_true(is.finite(b_est) && is.finite(c_est))
})
