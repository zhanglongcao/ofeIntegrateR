test_that("cv_krige_surface returns sane diagnostics", {
  sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 25,
                            seed = 4)
  cv <- cv_krige_surface(sim$point_samples, value = "point_obs")

  expect_named(cv, c("rmse", "r2", "n", "residuals"))
  expect_equal(cv$n, 25L)
  expect_length(cv$residuals, 25L)
  expect_gt(cv$rmse, 0)
  expect_lte(cv$r2, 1)
  expect_false(is.null(attr(cv, "variogram")))
})

test_that("r2 is negative when the surface carries no spatial signal", {
  # Pure noise at random locations: kriging cannot beat the mean
  set.seed(9)
  pts <- data.frame(col = runif(30, 1, 40), row = runif(30, 1, 40),
                    noise = rnorm(30))
  cv <- cv_krige_surface(pts, value = "noise")
  expect_lt(cv$r2, 0.2)
})

test_that("cv_krige_surface validates its inputs", {
  sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 25,
                            seed = 4)
  expect_error(cv_krige_surface(sim$point_samples, value = "nope"), "not found")
  expect_error(cv_krige_surface(sim$point_samples, value = "point_obs",
                                coords = c("lon", "lat")), "coords")
  expect_error(cv_krige_surface(sim$point_samples[1:3, ], value = "point_obs"),
               "[Aa]t least 4")
})
