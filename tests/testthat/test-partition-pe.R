planted <- function(seed = 7, step = 1.5, cut_at = 22) {
  sim <- simulate_ofe_trial(n_row = 40, n_col = 12, point_range = 4, seed = seed)
  d <- sim$grid
  d$dense_response <- d$dense_response + ifelse(d$row > cut_at, step, 0)
  d
}

test_that("a planted step is found, at about the right place", {
  d <- planted()
  z <- partition_pseudo_env(d, response = "dense_response", along = "row",
                            treat = "treat")
  p <- attr(z, "partition")
  expect_equal(p$n_zones, 2L)
  # The break should land near row 22, not at an arbitrary place.
  expect_lt(abs(p$breaks[2] - 22), 5)
  expect_true(is.factor(z$zone))
  expect_equal(nlevels(z$zone), 2L)
  expect_equal(nrow(z), nrow(d))
})

test_that("the derivation is reported, not just the answer", {
  z <- partition_pseudo_env(planted(), response = "dense_response",
                            along = "row", treat = "treat")
  p <- attr(z, "partition")
  expect_named(p, c("method", "along", "response", "n_zones", "range",
                    "min_zone_width", "spacing", "extent",
                    "max_zones_possible", "breaks", "zones", "bic", "profile"))
  expect_gt(p$range, 0)
  expect_gte(p$min_zone_width, p$spacing)
  expect_equal(sum(p$zones$n), nrow(z))
  expect_equal(nrow(p$profile), 40L)
  expect_true(all(diff(p$bic$n_zones) == 1))
})

test_that("the fitted range caps how finely the trial can be cut", {
  d <- planted()
  # Insisting on a zone narrower than the correlation range is refused rather
  # than granted, which is the whole point of estimating the range.
  z_wide <- partition_pseudo_env(d, response = "dense_response", along = "row",
                                 treat = "treat", min_zone_width = 30)
  expect_equal(attr(z_wide, "partition")$max_zones_possible, 1L)
  expect_equal(attr(z_wide, "partition")$n_zones, 1L)

  z_fine <- partition_pseudo_env(d, response = "dense_response", along = "row",
                                 treat = "treat", min_zone_width = 2,
                                 max_zones = 6)
  expect_gte(attr(z_fine, "partition")$max_zones_possible, 6L)
})

test_that("a fixed number of zones can be demanded, equal or fitted", {
  d <- planted()
  z3 <- partition_pseudo_env(d, response = "dense_response", along = "row",
                             treat = "treat", n_zones = 3,
                             min_zone_width = 5)
  expect_equal(nlevels(z3$zone), 3L)

  ze <- partition_pseudo_env(d, response = "dense_response", along = "row",
                             method = "equal", n_zones = 4)
  expect_equal(nlevels(ze$zone), 4L)
  expect_equal(unname(diff(attr(ze, "partition")$breaks)), rep(10, 4))
  expect_error(partition_pseudo_env(d, response = "dense_response",
                                    along = "row", method = "equal"),
               "needs `n_zones`")
  expect_error(partition_pseudo_env(d, response = "dense_response",
                                    along = "row", treat = "treat",
                                    n_zones = 8, min_zone_width = 20),
               "not reachable")
})

test_that("bad input is refused", {
  d <- planted()
  expect_error(partition_pseudo_env(d, response = "nope", along = "row"),
               "not found in `data`")
  expect_error(partition_pseudo_env(d, response = "dense_response",
                                    along = "treat"), "numeric")
})

test_that("adaptive_residual writes a formula both engines accept", {
  z <- partition_pseudo_env(planted(), response = "dense_response",
                            along = "row", treat = "treat")
  r <- adaptive_residual(z)
  expect_s3_class(r, "formula")
  expect_match(deparse(r), "dsum")
  expect_match(deparse(r), "ar1\\(row\\):ar1\\(col\\)")
  expect_equal(attr(r, "n_degenerate"), 0L)
  expect_equal(nrow(attr(r, "geometry")), 2L)
  # It has to be usable, not merely well formed.
  fit <- fit_ofe(dense_response ~ zone + zone:treat, data = z, residual = r)
  expect_s3_class(fit, "ofe_fit")
  expect_equal(length(fit$sections), 2L)
})

test_that("adaptive_residual demotes a zone that has no extent to id()", {
  d <- data.frame(row = rep(1:6, each = 4), col = rep(1:4, 6))
  d$zone <- factor(ifelse(d$row <= 5, "wide", "thin"))
  r <- adaptive_residual(d)
  expect_equal(attr(r, "n_degenerate"), 1L)
  geo <- attr(r, "geometry")
  expect_match(geo$struct[geo$zone == "thin"], "id\\(row\\)")
  expect_match(geo$struct[geo$zone == "wide"], "ar1\\(row\\)")
  expect_error(adaptive_residual(d, zone = "nope"), "not found in `data`")
})

test_that("zones derived from an independent layer also work", {
  # The documented way to avoid selecting zones on the response being analysed.
  sim <- simulate_ofe_trial(n_row = 40, n_col = 12, point_range = 4, seed = 2)
  d <- sim$grid
  d$em38 <- d$point_true + ifelse(d$row > 20, 2, 0)
  z <- partition_pseudo_env(d, response = "em38", along = "row")
  expect_equal(attr(z, "partition")$n_zones, 2L)
  expect_lt(abs(attr(z, "partition")$breaks[2] - 20), 4)
})
