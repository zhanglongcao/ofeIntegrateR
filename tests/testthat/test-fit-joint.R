test_that("fit_integrated_joint requires the asreml package", {
  skip_if(requireNamespace("asreml", quietly = TRUE),
          "asreml is installed; informative-error path not applicable")
  sim <- simulate_ofe_trial(n_row = 10, n_col = 8, n_point_samples = 6, seed = 6)
  expect_error(
    fit_integrated_joint(sim$grid, sim$point_samples,
                          response_dense = "dense_response",
                          response_point = "point_obs"),
    "asreml"
  )
})

test_that("fit_integrated_joint warns when too few points for cor_structure = 'unstructured'", {
  skip_if_not_installed("asreml")
  sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 6, seed = 7)
  expect_warning(
    fit_integrated_joint(sim$grid, sim$point_samples,
                          response_dense = "dense_response",
                          response_point = "point_obs",
                          cor_structure = "unstructured",
                          min_point_n_for_us = 30),
    "fewer than"
  )
})

test_that("fit_integrated_joint (cor_structure = 'independent') returns sensible treatment contrasts", {
  skip_if_not_installed("asreml")
  sim <- simulate_ofe_trial(n_row = 30, n_col = 21, n_treat = 3,
                             treat_effects = c(0, 0.8, 1.6),
                             n_point_samples = 25, seed = 11)
  fit <- fit_integrated_joint(sim$grid, sim$point_samples,
                               response_dense = "dense_response",
                               response_point = "point_obs")
  expect_s3_class(fit, "asreml")

  fx <- extract_fixed_effects(fit)
  b_est <- fx$estimate[grepl("treat_B$", fx$term)]
  c_est <- fx$estimate[grepl("treat_C$", fx$term)]
  expect_length(b_est, 1)
  expect_length(c_est, 1)
  expect_true(is.finite(b_est) && is.finite(c_est))
  expect_gt(c_est, b_est)
})
