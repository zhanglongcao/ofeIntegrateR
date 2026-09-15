test_that("denser sampling is required for shorter-range variation", {
  n <- sapply(c(120, 60, 30), function(r) {
    kriging_sample_interval(nugget = 0.2, psill = 0.8, range = r,
                            target_kse = 0.6, area_ha = 20)$n_samples
  })
  # Shorter range -> more samples, strictly
  expect_true(all(diff(n) > 0))
})

test_that("the interval scales with the variogram range", {
  a <- kriging_sample_interval(nugget = 0.2, psill = 0.8, range = 60,
                               target_kse = 0.6)
  b <- kriging_sample_interval(nugget = 0.2, psill = 0.8, range = 120,
                               target_kse = 0.6)
  expect_gt(b$interval, a$interval)
  # Doubling the range should roughly double the interval
  expect_equal(b$interval / a$interval, 2, tolerance = 0.15)
})

test_that("the nugget sets a floor no sampling density can beat", {
  s <- kriging_sample_interval(nugget = 0.5, psill = 0.5, range = 60,
                               target_kse = 0.9)
  expect_equal(s$kse_floor, sqrt(0.5) / sqrt(1.0), tolerance = 1e-8)
  # Nothing in the curve may go below the floor
  expect_true(all(s$curve$rel_kse >= s$kse_floor - 1e-6))
})

test_that("an unattainable target warns and returns NA rather than guessing", {
  expect_warning(
    s <- kriging_sample_interval(nugget = 0.5, psill = 0.5, range = 60,
                                 target_kse = 0.4, area_ha = 20),
    "unattainable")
  expect_true(is.na(s$interval))
  expect_true(is.na(s$n_samples))
})

test_that("the achieved error meets the target and the curve is monotone", {
  s <- kriging_sample_interval(nugget = 0.1, psill = 0.9, range = 60,
                               target_kse = 0.6, area_ha = 10)
  expect_lte(s$rel_kse, s$target_kse)
  expect_true(is.finite(s$n_samples))
  # Wider spacing is never more precise
  expect_true(all(diff(s$curve$rel_kse) >= -1e-8))
})

test_that("kriging_sample_interval validates its inputs", {
  expect_error(kriging_sample_interval(0.2, 0.8, range = -1), "positive")
  expect_error(kriging_sample_interval(-1, 0.8, range = 60), "non-negative")
  expect_error(kriging_sample_interval(0.2, 0.8, range = 60, target_kse = 1.5),
               "between 0 and 1")
})
