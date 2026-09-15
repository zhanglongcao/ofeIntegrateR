test_that("simulate_yield_monitor produces an irregular, gappy cloud", {
  sim <- simulate_yield_monitor(n_point_samples = 20, seed = 2)

  expect_named(sim, c("cloud", "point_samples", "true_effects"))
  expect_true(all(c("x", "y", "treat", "point_true", "yield") %in%
                    names(sim$cloud)))
  expect_gt(nrow(sim$cloud), 500)

  # Points are not on a lattice: y values within a pass vary continuously
  expect_gt(length(unique(sim$cloud$y)), nrow(sim$cloud) / 2)

  # The clipped corner leaves the paddock non-rectangular
  far <- sim$cloud[sim$cloud$x > 0.9 * max(sim$cloud$x), ]
  expect_lt(max(far$y), max(sim$cloud$y))
})

test_that("the cloud grids into a complete lattice with genuine gaps", {
  sim <- simulate_yield_monitor(n_point_samples = 20, seed = 2)
  g <- grid_dense_layer(sim$cloud, response = "yield", treat = "treat",
                        cell_size = 9)
  expect_equal(nrow(g), max(g$row) * max(g$col))
  expect_true(any(g$n_obs == 0))
})

test_that("strip_ends sampling leaves the middle of the trial unsampled", {
  se <- simulate_yield_monitor(n_point_samples = 12, point_design = "strip_ends",
                               seed = 2)
  rn <- simulate_yield_monitor(n_point_samples = 12, point_design = "random",
                               seed = 2)
  mid <- function(p, cloud) {
    lo <- stats::quantile(cloud$x, 0.3); hi <- stats::quantile(cloud$x, 0.7)
    mean(p$x > lo & p$x < hi)
  }
  expect_lt(mid(se$point_samples, se$cloud), mid(rn$point_samples, rn$cloud))
})

test_that("point samples sit at arbitrary coordinates, not cell centres", {
  sim <- simulate_yield_monitor(n_point_samples = 20, seed = 3)
  expect_true(all(c("x", "y", "point_obs") %in% names(sim$point_samples)))
  expect_false(all(sim$point_samples$x %% 1 == 0))
})
