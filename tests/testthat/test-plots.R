# Plot code is tested for running cleanly and for the decisions it makes about
# scale and colour, which are the parts that can be wrong without erroring.
skip_plot_setup <- function() {
  grDevices::pdf(NULL)
  withr_teardown <- grDevices::dev.cur()
  withr_teardown
}

test_that("the palettes are the validated ones and the right length", {
  expect_equal(length(ofe_palette("treatment", 4)), 4L)
  expect_equal(ofe_palette("treatment", 3),
               c("#2a78d6", "#1baf7a", "#eda100"))
  expect_equal(length(ofe_palette("zone", 6)), 6L)
  expect_equal(length(ofe_palette("surface", 32)), 32L)
  expect_equal(length(ofe_palette("residual", 9)), 9L)
  for (w in c("treatment", "zone", "surface", "residual")) {
    cols <- ofe_palette(w, 5)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)), info = w)
  }
  expect_named(ofe_palette("chrome"),
               c("surface", "ink", "ink2", "muted", "grid", "axis",
                 "midpoint"))
  # The categorical order is nested: asking for more never repaints the first.
  expect_equal(ofe_palette("treatment", 2), ofe_palette("treatment", 5)[1:2])
})

test_that("the diverging scale meets at a neutral, with opposite poles", {
  d <- ofe_palette("residual", 5)
  rgbv <- grDevices::col2rgb(d)
  # The midpoint must be near-grey: a hue there would invent a category.
  mid <- rgbv[, 3]
  expect_lt(diff(range(mid)), 12)
  # The poles are a cool and a warm colour, so the sign reads without a legend.
  expect_gt(rgbv["blue", 1], rgbv["red", 1])
  expect_gt(rgbv["red", 5], rgbv["blue", 5])
})

test_that("auto never picks a diverging scale for a merely negative variable", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- expand.grid(x = 1:6, y = 1:6)
  d$v <- seq(-5, 5, length.out = 36)          # crosses zero, but not signed
  key <- ofe_map(d, "v", legend = FALSE)
  # A diverging key would be symmetric about zero; a sequential one spans the
  # data. This is the check that the auto rule did not silently change back.
  expect_equal(key$zlim, range(d$v))
  key_div <- ofe_map(d, "v", type = "diverging", legend = FALSE)
  expect_equal(key_div$zlim, c(-5, 5))
})

test_that("ofe_map draws each scale type and reports its key", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  sim <- simulate_ofe_trial(n_row = 12, n_col = 8, seed = 1)
  seq_key <- ofe_map(sim$grid, "dense_response")
  expect_equal(length(seq_key$colours), 64L)
  expect_equal(length(seq_key$fill), nrow(sim$grid))

  cat_key <- ofe_map(sim$grid, "treat")
  expect_equal(cat_key$levels, levels(sim$grid$treat))
  expect_equal(length(cat_key$colours), nlevels(sim$grid$treat))

  expect_error(ofe_map(sim$grid, "nope"), "not found in `data`")
  expect_error(ofe_map(data.frame(v = 1:4), "v"), "Coordinates not found")
})

test_that("missing values are drawn as missing, not as zero", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- expand.grid(x = 1:5, y = 1:5)
  d$v <- 1
  d$v[3] <- NA
  key <- ofe_map(d, "v", na_col = "#123456", legend = FALSE)
  expect_equal(key$fill[3], "#123456")
  expect_false(any(key$fill[-3] == "#123456"))
})

test_that("the overlay runs inside the map's coordinate system", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  sim <- simulate_ofe_trial(n_row = 10, n_col = 6, seed = 1)
  seen <- NULL
  ofe_map(sim$grid, "dense_response",
          overlay = function() seen <<- graphics::par("usr"))
  expect_false(is.null(seen))
  # The overlay must see the map's limits, not the device default (0, 1).
  expect_gt(seen[2], 5)
})

test_that("label ink is whichever actually contrasts with the fill", {
  chrome <- ofe_palette("chrome")
  expect_equal(ofeIntegrateR:::.ofe_ink_on("#0d366b"),
               unname(chrome[["surface"]]))      # deep blue -> light ink
  expect_equal(ofeIntegrateR:::.ofe_ink_on("#eda100"),
               unname(chrome[["ink"]]))          # mid yellow -> dark ink
  expect_equal(ofeIntegrateR:::.ofe_ink_on("#ffffff"), unname(chrome[["ink"]]))
  expect_equal(ofeIntegrateR:::.ofe_ink_on("#000000"),
               unname(chrome[["surface"]]))
  expect_length(ofeIntegrateR:::.ofe_ink_on(c("#000000", "#ffffff")), 2L)
  expect_null(names(ofeIntegrateR:::.ofe_ink_on("#000000")))
})

test_that("plot methods run and return their object invisibly", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  d <- make_trial_design(c("N0", "N60", "N120"), n_rep = 3, seed = 1)
  expect_s3_class(d, "ofe_design")
  expect_invisible(plot(d))

  g <- expand.grid(row = 1:20, col = 1:14)
  g$x <- g$col * 10
  g$y <- g$row * 10
  set.seed(1)
  g$elevation <- 100 + 5 * (g$y > 100) + stats::rnorm(nrow(g), 0, 0.3)
  z <- partition_paddock(g, covariates = "elevation")
  expect_s3_class(z, "ofe_zones")
  expect_invisible(plot(z))

  sim <- simulate_ofe_trial(n_row = 30, n_col = 12, point_range = 5,
                            n_point_samples = 60, seed = 4)
  zz <- partition_pseudo_env(sim$grid, response = "dense_response",
                             along = "row", treat = "treat")
  expect_s3_class(zz, "ofe_zones")
  expect_invisible(plot(zz))

  v <- ofe_variogram(within(sim$point_samples, {x <- col; y <- row}),
                     value = "point_obs")
  expect_invisible(plot(v))
  expect_invisible(plot(v, show_all = TRUE))
})

test_that("fit diagnostics draw, and report the residual variogram back", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  sim <- simulate_ofe_trial(n_row = 20, n_col = 12, seed = 1)
  fit <- fit_ofe(dense_response ~ treat, data = sim$grid)
  out <- plot(fit)
  expect_s3_class(out$variogram, "ofe_variogram")
  # Panels can be chosen, and panel 4 is the only one that fits a variogram.
  expect_null(plot(fit, which = 2:3)$variogram)
  expect_null(plot(fit, which = integer(0)))
})

test_that("classing the pipeline outputs left them usable as data frames", {
  d <- make_trial_design(c("A", "B"), n_rep = 2, seed = 1)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d[d$treat == "A" & !is.na(d$treat), ]),
               sum(d$treat == "A", na.rm = TRUE))
  expect_true(is.data.frame(as.data.frame(d)))
})
