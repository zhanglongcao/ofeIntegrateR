# A field with known variogram parameters, sampled densely enough that the
# fit has something to work with.
sampled_field <- function(seed = 2, n = 150, range = 2) {
  sim <- simulate_ofe_trial(n_row = 30, n_col = 20, n_point_samples = n,
                            point_range = range, seed = seed)
  p <- sim$point_samples
  p$x <- p$col
  p$y <- p$row
  p
}

test_that("a variogram is fitted and its parts hang together", {
  v <- ofe_variogram(sampled_field(), value = "point_obs")
  expect_s3_class(v, "ofe_variogram")
  expect_equal(v$sill, v$nugget + v$psill)
  expect_equal(v$nugget_ratio, v$nugget / v$sill)
  expect_gte(v$nugget, 0)
  expect_gt(v$psill, 0)
  expect_gt(v$practical_range, 0)
  expect_true(all(c("h", "gamma", "n_pair") %in% names(v$empirical)))
  expect_true(all(v$empirical$n_pair >= 2))
  expect_true(all(v$empirical$h <= v$cutoff))
  expect_equal(v$n, nrow(sampled_field()))
})

test_that("the fitted range tracks the range the field was simulated with", {
  short <- ofe_variogram(sampled_field(range = 2), value = "point_obs")
  long <- ofe_variogram(sampled_field(range = 8), value = "point_obs")
  expect_lt(short$practical_range, long$practical_range)
})

test_that("a range beyond the sampled lags is flagged, not quietly reported", {
  # A range long relative to the paddock: the variogram cannot reach its sill
  # inside a cutoff of a third of the largest distance.
  v <- ofe_variogram(sampled_field(seed = 1, range = 6, n = 60), value = "point_obs")
  expect_false(v$range_identified)
  expect_gt(v$practical_range, v$max_lag)
  expect_output(print(v), "extrapolation")

  ok <- ofe_variogram(sampled_field(range = 2, n = 150), value = "point_obs")
  expect_true(ok$range_identified)
  expect_lte(ok$practical_range, ok$max_lag)
})

test_that("print classifies spatial dependence the usual way", {
  v <- ofe_variogram(sampled_field(), value = "point_obs")
  out <- utils::capture.output(print(v))
  cls <- if (v$nugget_ratio < 0.25) "strong" else
    if (v$nugget_ratio <= 0.75) "moderate" else "weak"
  expect_true(any(grepl(cls, out, fixed = TRUE)))
})

test_that("a trend is removed before the variogram, not left to inflate it", {
  p <- sampled_field()
  p$treat <- factor(rep(c("A", "B"), length.out = nrow(p)))
  p$point_obs <- p$point_obs + ifelse(p$treat == "B", 6, 0)
  raw <- ofe_variogram(p, value = "point_obs")
  adj <- ofe_variogram(p, value = "point_obs", trend = ~ treat)
  # The treatment step is pure noise as far as space is concerned, so leaving
  # it in inflates the sill; removing it must bring the sill back down.
  expect_gt(raw$sill, adj$sill * 2)
})

test_that("several models can be compared and the best is kept", {
  v <- ofe_variogram(sampled_field(), value = "point_obs",
                     model = c("exponential", "spherical", "gaussian"))
  expect_equal(length(v$fits), 3L)
  expect_equal(v$sse, min(vapply(v$fits, function(z) z$sse, numeric(1))))
  expect_true(v$model %in% c("exponential", "spherical", "gaussian"))
  one <- ofe_variogram(sampled_field(), value = "point_obs",
                       model = "spherical")
  expect_equal(one$model, "spherical")
  expect_gte(one$sse, v$sse)
})

test_that("variograms refuse input they cannot describe", {
  p <- sampled_field()
  expect_error(ofe_variogram(p, value = "nope"), "not found in `data`")
  expect_error(ofe_variogram(p[1:8, ], value = "point_obs"), "too few")
  expect_error(ofe_variogram(data.frame(z = stats::rnorm(50)), value = "z"),
               "Coordinates not found")
})

test_that("kriging_sample_interval takes a fitted variogram directly", {
  v <- ofe_variogram(sampled_field(range = 2, n = 150), value = "point_obs")
  from_obj <- kriging_sample_interval(v, target_kse = 0.6, area_ha = 4)
  by_hand <- kriging_sample_interval(nugget = v$nugget, psill = v$psill,
                                     range = v$practical_range, model = "Sph",
                                     target_kse = 0.6, area_ha = 4)
  if (v$model == "spherical") expect_equal(from_obj, by_hand)
  expect_true(is.finite(from_obj$kse_floor))
  expect_error(kriging_sample_interval(v, psill = 1), "not both")
  expect_error(kriging_sample_interval(nugget = 0.2), "Supply `psill`")
})

test_that("an unidentified range warns when it is turned into a sample plan", {
  v <- ofe_variogram(sampled_field(seed = 1, range = 6, n = 60), value = "point_obs")
  expect_warning(kriging_sample_interval(v, target_kse = 0.7),
                 "extrapolation")
})
