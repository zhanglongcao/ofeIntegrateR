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

test_that("fit_integrated_kriged (engine = 'gls') returns sensible treatment contrasts", {
  skip_if_not_installed("nlme")
  sim <- simulate_ofe_trial(n_row = 30, n_col = 21, n_treat = 3,
                             treat_effects = c(0, 0.8, 1.6),
                             n_point_samples = 25, seed = 11)
  krieged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")

  fit <- fit_integrated_kriged(krieged, response = "dense_response",
                                treat = "treat", covariate = "point_obs_kriged",
                                row = "row", col = "col", engine = "gls")
  expect_s3_class(fit, "gls")

  fx <- extract_fixed_effects(fit)
  b_est <- fx$estimate[fx$term == "treatB"]
  c_est <- fx$estimate[fx$term == "treatC"]
  expect_true(is.finite(b_est) && is.finite(c_est))
  expect_gt(b_est, 0)
  expect_gt(c_est, b_est)
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

test_that("several covariates all enter the model (none silently dropped)", {
  sim <- simulate_ofe_trial(n_row = 20, n_col = 9, n_point_samples = 15, seed = 42)
  g <- sim$grid
  set.seed(42)
  g$soil_n <- stats::rnorm(nrow(g))
  g$soil_p <- stats::rnorm(nrow(g))

  fit1 <- fit_integrated_kriged(g, "dense_response", "treat", "soil_n",
                                engine = "lm")
  fit2 <- fit_integrated_kriged(g, "dense_response", "treat",
                                c("soil_n", "soil_p"), engine = "lm")

  expect_true("soil_n" %in% names(stats::coef(fit1)))
  expect_false("soil_p" %in% names(stats::coef(fit1)))

  # Both covariates must be present when both are requested.
  expect_true(all(c("soil_n", "soil_p") %in% names(stats::coef(fit2))))
  expect_equal(length(attr(stats::terms(fit2), "term.labels")), 3L)
})

test_that("gls engine tolerates missing cells, as gridded trial data contain", {
  sim <- simulate_ofe_trial(n_row = 14, n_col = 9, n_point_samples = 12, seed = 7)
  kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")

  # Empty cells and treatment-boundary cells arrive as NA from grid_dense_layer()
  kr$dense_response[c(3, 17, 42)] <- NA_real_

  expect_error(
    fit <- fit_integrated_kriged(kr, "dense_response", "treat",
                                 "point_obs_kriged", engine = "gls"),
    NA)
  fx <- extract_fixed_effects(fit)
  expect_true(any(grepl("treatC$", fx$term)))
  expect_true(all(is.finite(fx$estimate)))
})

test_that("covariate = NULL fits the dense-layer-only baseline on every engine", {
  sim <- simulate_ofe_trial(n_row = 14, n_col = 9, n_point_samples = 12, seed = 3)
  kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")

  for (eng in c("lm", "gls")) {
    fit <- fit_integrated_kriged(kr, "dense_response", "treat", covariate = NULL,
                                 engine = eng)
    fx <- extract_fixed_effects(fit)
    expect_true(any(grepl("treatC$", fx$term)), info = eng)
    # No covariate term should appear
    expect_false(any(grepl("kriged", fx$term)), info = eng)
    expect_true(all(is.finite(fx$estimate)), info = eng)
  }
})

test_that("the lme engine fits random effects with a spatial residual", {
  sim <- simulate_ofe_trial(n_row = 24, n_col = 12, n_point_samples = 20,
                            seed = 3)
  kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
  kr$rep <- factor(ceiling(kr$row / 8))

  fit <- fit_integrated_kriged(kr, "dense_response", "treat",
                               "point_obs_kriged", random = ~ rep,
                               engine = "lme")
  expect_s3_class(fit, "lme")

  fx <- extract_fixed_effects(fit)
  est <- fx$estimate[grepl("treatC$", fx$term)][1]
  expect_equal(est, 1.6, tolerance = 0.5)
  expect_true(all(is.finite(fx$se)))
})

test_that("the asreml spelling of `random` works on lme too", {
  sim <- simulate_ofe_trial(n_row = 24, n_col = 12, n_point_samples = 20,
                            seed = 3)
  kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
  kr$rep <- factor(ceiling(kr$row / 8))

  short <- fit_integrated_kriged(kr, "dense_response", "treat",
                                 "point_obs_kriged", random = ~ rep,
                                 engine = "lme")
  explicit <- fit_integrated_kriged(kr, "dense_response", "treat",
                                    "point_obs_kriged", random = ~ 1 | rep,
                                    engine = "lme")
  expect_equal(extract_fixed_effects(short)$estimate,
               extract_fixed_effects(explicit)$estimate, tolerance = 1e-8)
})

test_that("engines that cannot fit random effects say so instead of ignoring them", {
  sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 15,
                            seed = 5)
  kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
  kr$rep <- factor(ceiling(kr$row / 5))

  for (eng in c("gls", "lm")) {
    expect_error(
      fit_integrated_kriged(kr, "dense_response", "treat", "point_obs_kriged",
                            random = ~ rep, engine = eng),
      "cannot fit random effects", info = eng)
  }
})

test_that("lme requires a random formula and points at gls when there is none", {
  sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 15,
                            seed = 5)
  kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
  expect_error(
    fit_integrated_kriged(kr, "dense_response", "treat", "point_obs_kriged",
                          engine = "lme"),
    "needs a `random` formula")
})

test_that("ambiguous multi-factor random formulas are refused, not guessed", {
  sim <- simulate_ofe_trial(n_row = 24, n_col = 12, n_point_samples = 20,
                            seed = 3)
  kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
  kr$rep <- factor(ceiling(kr$row / 8))
  kr$blk <- factor(ceiling(kr$col / 4))
  expect_error(
    fit_integrated_kriged(kr, "dense_response", "treat", "point_obs_kriged",
                          random = ~ rep + blk, engine = "lme"),
    "nesting made explicit")
})
