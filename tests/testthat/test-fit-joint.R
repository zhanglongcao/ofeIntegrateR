test_that("fit_integrated_joint requires the asreml package", {
  skip_if(requireNamespace("asreml", quietly = TRUE),
          "asreml is installed; informative-error path not applicable")
  sim <- simulate_ofe_trial(n_row = 10, n_col = 8, n_point_samples = 12, seed = 6)
  expect_error(
    fit_integrated_joint(sim$grid, sim$point_samples,
                         response_dense = "dense_response",
                         response_point = "point_obs"),
    "asreml"
  )
})

test_that("the removed cor_structure argument gets an informative error", {
  sim <- simulate_ofe_trial(n_row = 10, n_col = 8, n_point_samples = 12, seed = 6)
  expect_error(
    fit_integrated_joint(sim$grid, sim$point_samples,
                         response_dense = "dense_response",
                         response_point = "point_obs",
                         cor_structure = "independent"),
    "removed"
  )
})

test_that("sommer refuses the spatial residual it cannot fit", {
  sim <- simulate_ofe_trial(n_row = 10, n_col = 8, n_point_samples = 12, seed = 6)
  expect_error(
    fit_integrated_joint(sim$grid, sim$point_samples,
                         response_dense = "dense_response",
                         response_point = "point_obs",
                         engine = "sommer", residual = "ar1"),
    "sommer cannot"
  )
})

test_that("point samples off the lattice are dropped with a warning", {
  sim <- simulate_ofe_trial(n_row = 10, n_col = 8, n_point_samples = 12, seed = 6)
  pts <- sim$point_samples
  pts$row[1] <- 999
  expect_warning(
    wide <- build_joint_wide(sim$grid, pts, "dense_response", "point_obs",
                             "row", "col"),
    "fall on no cell"
  )
  expect_equal(sum(!is.na(wide$point_obs)), 11)
})

test_that("several samples in one cell are averaged", {
  sim <- simulate_ofe_trial(n_row = 10, n_col = 8, n_point_samples = 4, seed = 6)
  pts <- rbind(sim$point_samples, sim$point_samples[1, ])
  pts$point_obs[5] <- pts$point_obs[1] + 2
  expect_message(
    wide <- build_joint_wide(sim$grid, pts, "dense_response", "point_obs",
                             "row", "col"),
    "mean"
  )
  hit <- wide$row == pts$row[1] & wide$col == pts$col[1]
  expect_equal(wide$point_obs[hit], pts$point_obs[1] + 1)
})

test_that("a sparse overlap is warned about, and an incomplete lattice refused", {
  sim <- simulate_ofe_trial(n_row = 10, n_col = 8, n_point_samples = 4, seed = 6)
  holed <- sim$grid[-nrow(sim$grid), ]
  pts <- sim$point_samples[paste(sim$point_samples$row, sim$point_samples$col) %in%
                             paste(holed$row, holed$col), ]
  expect_warning(
    expect_error(
      fit_integrated_joint(holed, pts,
                           response_dense = "dense_response",
                           response_point = "point_obs", min_overlap = 10),
      "complete row-by-column lattice|asreml"),
    "cell\\(s\\) carry both layers"
  )
})

test_that("asreml fits the spatial bivariate model and estimates the cross-covariance", {
  skip_if_not_installed("asreml")
  sim <- simulate_ofe_trial(n_row = 30, n_col = 21, n_treat = 3,
                            treat_effects = c(0, 0.8, 1.6),
                            n_point_samples = 40, seed = 11)
  fit <- fit_integrated_joint(sim$grid, sim$point_samples,
                              response_dense = "dense_response",
                              response_point = "point_obs")
  expect_s3_class(fit, "asreml")
  expect_true(fit$converge)

  vc <- summary(fit)$varcomp
  expect_true(any(grepl("trait_point_obs:dense_response", rownames(vc))))
  expect_true(any(grepl("!row!cor", rownames(vc))))

  fx <- extract_fixed_effects(fit)
  b_est <- fx$estimate[fx$term == "at(trait, dense_response):treat_B"]
  c_est <- fx$estimate[fx$term == "at(trait, dense_response):treat_C"]
  expect_length(b_est, 1)
  expect_length(c_est, 1)
  expect_gt(c_est, b_est)
})

test_that("the joint model is not the univariate model in disguise", {
  # The old diag(layer):id(unit) form reproduced a univariate fit exactly.
  # With a residual cross-covariance, the sparse layer's contrast must move
  # when the dense layer is added.
  skip_if_not_installed("asreml")
  sim <- simulate_ofe_trial(n_row = 30, n_col = 21, n_treat = 3,
                            treat_effects = c(0, 0.8, 1.6),
                            n_point_samples = 60, seed = 3)
  fit <- fit_integrated_joint(sim$grid, sim$point_samples,
                              response_dense = "dense_response",
                              response_point = "point_obs",
                              residual = "id")
  wide <- build_joint_wide(sim$grid, sim$point_samples, "dense_response",
                           "point_obs", "row", "col")
  uni <- asreml::asreml(point_obs ~ treat, data = wide[!is.na(wide$point_obs), ],
                        trace = FALSE)
  fx <- extract_fixed_effects(fit)
  fu <- extract_fixed_effects(uni)
  expect_false(isTRUE(all.equal(
    fx$estimate[fx$term == "at(trait, point_obs):treat_C"],
    fu$estimate[fu$term == "treat_C"])))
})

test_that("a block factor is fitted within each trait", {
  skip_if_not_installed("asreml")
  sim <- simulate_ofe_trial(n_row = 30, n_col = 21, n_point_samples = 40, seed = 11)
  g <- sim$grid
  g$rep <- factor((g$row - 1) %/% 10 + 1)
  fit <- fit_integrated_joint(g, sim$point_samples,
                              response_dense = "dense_response",
                              response_point = "point_obs", block = "rep")
  fx <- extract_fixed_effects(fit)
  expect_true(any(grepl("trait_dense_response:rep_2", fx$term)))
  expect_true(any(grepl("trait_point_obs:rep_2", fx$term)))
})

test_that("sommer refuses a sparse layer, which it would impute", {
  skip_if_not_installed("sommer")
  sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 30, seed = 11)
  expect_error(
    fit_integrated_joint(sim$grid, sim$point_samples,
                         response_dense = "dense_response",
                         response_point = "point_obs", engine = "sommer"),
    "both layers observed in every cell"
  )
})

test_that("sommer fits two complete layers and matches asreml", {
  skip_if_not_installed("sommer")
  sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_treat = 3,
                            treat_effects = c(0, 0.8, 1.6), seed = 11)
  both <- sim$grid[, c("row", "col", "point_true")]
  fit <- fit_integrated_joint(sim$grid, both,
                              response_dense = "dense_response",
                              response_point = "point_true",
                              engine = "sommer", min_overlap = 0)
  expect_s3_class(fit, "mmer")
  fx <- extract_fixed_effects(fit)
  c_s <- fx$estimate[fx$term == "dense_response:treatC"]
  expect_gt(c_s, fx$estimate[fx$term == "dense_response:treatB"])

  skip_if_not_installed("asreml")
  fa <- fit_integrated_joint(sim$grid, both,
                             response_dense = "dense_response",
                             response_point = "point_true",
                             residual = "id", min_overlap = 0)
  fxa <- extract_fixed_effects(fa)
  expect_equal(c_s, fxa$estimate[fxa$term == "at(trait, dense_response):treat_C"],
               tolerance = 1e-3)
})
