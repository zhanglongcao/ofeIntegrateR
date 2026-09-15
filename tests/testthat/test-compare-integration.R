test_that("compare_integration reports both models on the same terms", {
  sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 20,
                            seed = 8)
  kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
  cmp <- compare_integration(kr, response = "dense_response", treat = "treat",
                             covariate = "point_obs_kriged", engine = "lm")

  expect_true(all(c("term", "baseline_estimate", "integrated_estimate",
                    "estimate_change", "se_change") %in% names(cmp)))
  expect_true(any(grepl("treatC$", cmp$term)))

  # The covariate appears only in the integrated fit
  cov_row <- cmp[cmp$term == "point_obs_kriged", ]
  expect_equal(nrow(cov_row), 1L)
  expect_true(is.na(cov_row$baseline_estimate))
  expect_true(is.finite(cov_row$integrated_estimate))

  # Both fits are returned for inspection
  m <- attr(cmp, "models")
  expect_named(m, c("baseline", "integrated"))
})

test_that("compare_integration recovers the treatment effect", {
  sim <- simulate_ofe_trial(n_row = 30, n_col = 12, treat_effects = c(0, 0.8, 1.6),
                            n_point_samples = 30, seed = 12)
  kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
  cmp <- compare_integration(kr, "dense_response", "treat", "point_obs_kriged",
                             engine = "lm")
  est <- cmp$integrated_estimate[grepl("treatC$", cmp$term)]
  expect_equal(est, 1.6, tolerance = 0.3)
})

test_that("compare_integration requires something to compare", {
  sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 20,
                            seed = 8)
  expect_error(compare_integration(sim$grid, "dense_response", "treat",
                                   covariate = NULL, engine = "lm"),
               "at least one column")
})
