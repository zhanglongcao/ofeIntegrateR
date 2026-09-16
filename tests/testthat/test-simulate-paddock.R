test_that("a paddock comes back with its layers, zones and a trial", {
  p <- simulate_paddock(n_row = 24, n_col = 18, n_zones = 2,
                        treatments = c("A", "B", "C"), n_rep = 3, seed = 1)
  expect_true(all(c("row", "col", "x", "y", "yield_potential", "elevation",
                    "clay", "zone", "treat", "rep", "plot", "yield") %in%
                    names(p)))
  expect_equal(nrow(p), 24L * 18L)
  expect_s3_class(p, "ofe_zones")
  expect_equal(nlevels(p$zone), 2L)
  expect_false(anyNA(p$yield_potential))
})

test_that("covariates are correlated with yield potential to order", {
  # The point of the covariates: unrelated ones are not worth zoning on and
  # identical ones are an unrealistically easy test.
  p <- simulate_paddock(n_row = 40, n_col = 30, covariate_cor = 0.8, seed = 2)
  r <- attr(p, "truth")$covariate_cor
  expect_equal(unname(r[["elevation"]]), 0.8, tolerance = 0.15)
  expect_equal(unname(r[["clay"]]), 0.8, tolerance = 0.15)

  indep <- simulate_paddock(n_row = 40, n_col = 30, covariate_cor = 0, seed = 2)
  expect_lt(abs(attr(indep, "truth")$covariate_cor[["elevation"]]), 0.25)

  # Per-layer settings override the default, including the correlation.
  custom <- simulate_paddock(n_row = 30, n_col = 20, seed = 3,
                             covariates = list(em38 = list(mean = 50, sd = 10,
                                                           cor = 0.9)))
  expect_true("em38" %in% names(custom))
  expect_false("elevation" %in% names(custom))
  expect_equal(mean(custom$em38), 50, tolerance = 0.5)
  expect_equal(stats::sd(custom$em38), 10, tolerance = 0.5)
  expect_gt(attr(custom, "truth")$covariate_cor[["em38"]], 0.7)
})

test_that("zone effects scale the treatment response, zone by zone", {
  p <- simulate_paddock(n_row = 30, n_col = 24, n_zones = 2,
                        treatments = c("N0", "N120"),
                        treat_effect = c(0, 1), zone_effect = c(1, 3),
                        noise_sd = 0.05, seed = 4)
  gap <- tapply(p$yield, list(p$zone, p$treat), mean, na.rm = TRUE)
  d1 <- gap["1", "N120"] - gap["1", "N0"]
  d2 <- gap["2", "N120"] - gap["2", "N0"]
  expect_equal(unname(d1), 1, tolerance = 0.35)
  expect_equal(unname(d2 / d1), 3, tolerance = 0.9)

  # One zone means no interaction is there to be found, which is the case that
  # most needs to be testable.
  flat <- simulate_paddock(n_row = 24, n_col = 18, treatments = c("A", "B"),
                           seed = 5)
  expect_equal(attr(flat, "truth")$n_zones, 1L)
  expect_equal(nlevels(flat$zone), 1L)
})

test_that("unequal zones are allowed, and the breaks are recorded", {
  p <- simulate_paddock(n_row = 40, n_col = 18, n_zones = 3,
                        zone_props = c(0.5, 0.3, 0.2), seed = 6)
  share <- as.numeric(table(p$zone)) / nrow(p)
  expect_equal(share, c(0.5, 0.3, 0.2), tolerance = 0.05)
  expect_equal(length(attr(p, "truth")$zone_breaks), 4L)
  expect_error(simulate_paddock(n_zones = 2, zone_props = c(1, 1, 1)),
               "one element per zone")
})

test_that("the trial is balanced across the zones it is meant to be analysed in", {
  p <- simulate_paddock(n_row = 30, n_col = 24, n_zones = 2,
                        treatments = c("A", "B", "C"), n_rep = 4, seed = 7)
  tab <- table(p$zone, p$treat)
  # Strips run the length of the paddock and zones cut across them, so every
  # treatment appears in every zone -- the property a PE analysis needs.
  expect_true(all(tab > 0))
  expect_equal(length(unique(as.vector(tab))), 1L)
})

test_that("a paddock without treatments is still a paddock", {
  p <- simulate_paddock(n_row = 20, n_col = 14, seed = 8)
  expect_false("yield" %in% names(p))
  expect_false("treat" %in% names(p))
  expect_true("yield_potential" %in% names(p))
  expect_null(attr(p, "truth")$treat_effect)
})

test_that("the paddock is reproducible and refuses impossible trials", {
  a <- simulate_paddock(n_row = 20, n_col = 14, treatments = c("A", "B"),
                        seed = 11)
  b <- simulate_paddock(n_row = 20, n_col = 14, treatments = c("A", "B"),
                        seed = 11)
  expect_equal(a$yield, b$yield)
  expect_error(simulate_paddock(n_row = 20, n_col = 4,
                                treatments = c("A", "B", "C"), n_rep = 4),
               "too narrow")
  expect_error(simulate_paddock(treatments = c("A", "B"),
                                treat_effect = 1:3), "per treatment")
  expect_error(simulate_paddock(n_zones = 2, treatments = c("A", "B"),
                                zone_effect = 1), "per zone")
})

test_that("the zoning functions can be scored against the simulated truth", {
  p <- simulate_paddock(n_row = 40, n_col = 24, n_zones = 2,
                        covariate_cor = 0.9, yield_range = 40, seed = 12)
  z <- partition_paddock(p, covariates = c("elevation", "clay"), k = 2)
  # Not a claim that it recovers the bands -- the covariates are correlated
  # with yield potential, not with the zone labels, which are a separate
  # structure. What matters is that the comparison is now possible at all.
  expect_equal(dim(table(z$zone, p$zone)), c(2L, 2L))
  expect_s3_class(p, "ofe_zones")
})
