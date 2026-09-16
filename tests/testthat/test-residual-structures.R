test_that("structure formulae parse into the components they name", {
  d <- data.frame(row = rep(1:4, each = 3), col = rep(1:3, 4),
                  zone = factor(rep(c("a", "b"), each = 6)))
  spec <- ofeIntegrateR:::.build_residual_spec(~ ar1(row):ar1(col), d)
  expect_equal(length(spec$sections), 1L)
  expect_equal(spec$ptab$type, c("cor", "cor"))
  expect_equal(spec$ptab$name, c("R!row!cor", "R!col!cor"))

  # A bare name is id(), as in asreml.
  expect_equal(ofeIntegrateR:::.build_residual_spec(~ row:col, d)$n_par, 0L)

  ex <- ofeIntegrateR:::.build_residual_spec(~ exp(row):id(col), d)
  expect_equal(ex$ptab$type, "range")

  dg <- ofeIntegrateR:::.build_residual_spec(~ diag(zone):ar1(row), d)
  expect_equal(dg$ptab$type, c("var", "cor"))
})

test_that("dsum() splits the data into independent sections", {
  d <- data.frame(row = rep(1:4, each = 3), col = rep(1:3, 4),
                  zone = factor(rep(c("a", "b"), each = 6)))
  spec <- ofeIntegrateR:::.build_residual_spec(
    ~ dsum(~ ar1(row):ar1(col) | zone), d)
  expect_equal(length(spec$sections), 2L)
  expect_equal(spec$sections[[1]]$idx, 1:6)
  expect_equal(spec$sections[[2]]$idx, 7:12)
  # Two correlations per section, plus one relative variance for the second.
  expect_equal(spec$n_par, 5L)
  expect_true("zone!b!var" %in% spec$ptab$name)

  # Two dsum() terms with different structures, as adaptive_residual() writes.
  spec2 <- ofeIntegrateR:::.build_residual_spec(
    ~ dsum(~ ar1(row):ar1(col) | zone, levels = c("a")) +
      dsum(~ ar1(row):id(col) | zone, levels = c("b")), d)
  expect_equal(length(spec2$sections), 2L)
  expect_equal(sum(spec2$ptab$type == "cor"), 3L)
})

test_that("unusable structures are refused rather than fitted quietly", {
  d <- data.frame(row = rep(1:4, each = 3), col = rep(1:3, 4),
                  zone = factor(rep(c("a", "b"), each = 6)))
  expect_error(ofeIntegrateR:::.build_residual_spec(~ ar2(row), d),
               "Unsupported residual structure")
  expect_error(ofeIntegrateR:::.build_residual_spec(~ ar1(nope), d),
               "not found in `data`")
  expect_error(ofeIntegrateR:::.build_residual_spec(
    ~ dsum(~ ar1(row) | zone, levels = c("z")), d), "not present")
  expect_error(ofeIntegrateR:::.build_residual_spec(
    ~ dsum(~ ar1(row):ar1(col) | zone, levels = c("a")), d),
    "in no residual section")
  expect_error(ofeIntegrateR:::.build_residual_spec(~ ar1(row) + ar1(col), d),
               "sum of `dsum\\(\\)` terms")
})

test_that("an AR1 dimension with one level falls back to id() with a warning", {
  d <- data.frame(row = rep(1:6, each = 1), col = 1L,
                  zone = factor(rep(c("a", "b"), each = 3)))
  expect_warning(
    spec <- ofeIntegrateR:::.build_residual_spec(~ ar1(row):ar1(col), d),
    "one level")
  expect_equal(spec$n_par, 1L)
})

test_that("AR1 lags are counted on the whole trial, not within a section", {
  # Section b holds rows 5 and 7: the lag between them is 2, not 1, because a
  # row is missing. Getting this wrong quietly biases every AR1 estimate.
  d <- data.frame(row = c(1, 2, 3, 5, 7), col = 1L,
                  zone = factor(c("a", "a", "a", "b", "b")))
  spec <- ofeIntegrateR:::.build_residual_spec(
    ~ dsum(~ ar1(row):id(col) | zone), d)
  D <- spec$sections[[2]]$comps[[1]]$D
  expect_equal(D[1, 2], 2)
})
