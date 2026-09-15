# Build a small irregular "harvester pass" cloud over a 100 x 60 m trial.
make_cloud <- function(seed = 1) {
  set.seed(seed)
  pass_y <- rep(seq(2, 58, by = 4), each = 120)
  d <- data.frame(
    x = rep(seq(1, 99, length.out = 120), times = length(unique(pass_y))),
    y = pass_y + stats::rnorm(length(pass_y), 0, 0.3)
  )
  d$treat <- cut(d$x, breaks = c(0, 33, 66, 100), labels = c("A", "B", "C"))
  d$yield <- as.numeric(d$treat) + stats::rnorm(nrow(d), 0, 0.5)
  d
}

test_that("grid_dense_layer returns a complete rectangular lattice", {
  g <- grid_dense_layer(make_cloud(), response = "yield", cell_size = 5)

  expect_true(all(c("row", "col", "x_centre", "y_centre", "yield", "n_obs") %in% names(g)))
  # Complete lattice: every row x col combination present exactly once.
  expect_equal(nrow(g), max(g$row) * max(g$col))
  expect_equal(anyDuplicated(paste(g$row, g$col)), 0L)
  expect_equal(sort(unique(g$col)), seq_len(max(g$col)))
  expect_equal(sort(unique(g$row)), seq_len(max(g$row)))
})

test_that("empty cells are kept as NA so the lattice stays estimable", {
  cloud <- make_cloud()
  # Passes are 4 m apart, so 1 m cells leave many rows with no observations.
  g <- grid_dense_layer(cloud, response = "yield", cell_size = 1)

  expect_true(any(g$n_obs == 0L))
  expect_true(all(is.na(g$yield[g$n_obs == 0L])))
  expect_equal(nrow(g), max(g$row) * max(g$col))

  dropped <- grid_dense_layer(cloud, response = "yield", cell_size = 1,
                              keep_empty = FALSE)
  expect_true(nrow(dropped) < nrow(g))
  expect_true(all(dropped$n_obs > 0L))
})

test_that("cell aggregation and counts are correct", {
  cloud <- make_cloud()
  g <- grid_dense_layer(cloud, response = "yield", cell_size = 10)

  expect_equal(sum(g$n_obs), nrow(cloud))

  # Recompute one cell by hand.
  gr <- attr(g, "grid")
  ci <- floor((cloud$x - gr$x0) / gr$cell_size) + 1L
  ri <- floor((cloud$y - gr$y0) / gr$cell_size) + 1L
  hit <- ci == 2L & ri == 3L
  target <- g$yield[g$col == 2L & g$row == 3L]
  expect_equal(target, mean(cloud$yield[hit]), tolerance = 1e-10)
})

test_that("n_min blanks under-supported cells without deleting them", {
  g <- grid_dense_layer(make_cloud(), response = "yield", cell_size = 1,
                        n_min = 5L)
  expect_true(all(is.na(g$yield[g$n_obs < 5L])))
  expect_true(all(!is.na(g$yield[g$n_obs >= 5L])))
})

test_that("treatment is assigned by majority with a purity score", {
  g <- grid_dense_layer(make_cloud(), response = "yield", treat = "treat",
                        cell_size = 5)

  expect_true(is.factor(g$treat))
  expect_setequal(levels(g$treat), c("A", "B", "C"))
  expect_true(all(g$treat_purity[g$n_obs > 0] > 0 &
                    g$treat_purity[g$n_obs > 0] <= 1))
  # Cells wholly inside a strip are pure; boundary cells need not be.
  expect_true(mean(g$treat_purity == 1, na.rm = TRUE) > 0.7)
})

test_that("a non-default aggregation function is honoured", {
  cloud <- make_cloud()
  gm <- grid_dense_layer(cloud, response = "yield", cell_size = 10, fun = mean)
  gd <- grid_dense_layer(cloud, response = "yield", cell_size = 10, fun = stats::median)
  expect_false(isTRUE(all.equal(gm$yield, gd$yield)))
})

test_that("grid_dense_layer validates its inputs", {
  cloud <- make_cloud()
  expect_error(grid_dense_layer(cloud, response = "nope"), "not found")
  expect_error(grid_dense_layer(cloud, x = "lon", response = "yield"), "not found")
  expect_error(grid_dense_layer(cloud, response = "yield", cell_size = 0),
               "positive")
  bad <- cloud; bad$x <- NA_real_
  expect_error(grid_dense_layer(bad, response = "yield"), "finite coordinates")
})
