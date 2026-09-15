test_that("make_trial_grid lays out a complete, balanced lattice", {
  g <- make_trial_grid(width = 240, height = 108, cell_size = 9,
                       treatments = c("A", "B", "C"), strip_width = 18)

  expect_equal(nrow(g), max(g$row) * max(g$col))
  expect_setequal(names(g), c("row", "col", "x", "y", "treat"))
  expect_equal(max(g$col), 240 %/% 9)
  expect_equal(max(g$row), 108 %/% 9)
  # Six strips of 18 m over 108 m, cycling A B C twice: perfectly balanced
  expect_equal(length(unique(table(g$treat))), 1L)
})

test_that("strips run along the requested axis", {
  gx <- make_trial_grid(200, 100, cell_size = 10, strip_width = 20, along = "x")
  gy <- make_trial_grid(200, 100, cell_size = 10, strip_width = 20, along = "y")

  # Strips along x vary with y, and vice versa
  expect_equal(length(unique(gx$treat[gx$row == 1])), 1L)
  expect_equal(length(unique(gy$treat[gy$col == 1])), 1L)
  expect_gt(length(unique(gx$treat[gx$col == 1])), 1L)
  expect_gt(length(unique(gy$treat[gy$row == 1])), 1L)
})

test_that("the grid feeds place_point_samples directly", {
  g <- make_trial_grid(240, 108, cell_size = 9, strip_width = 18)
  set.seed(1)
  p <- place_point_samples(g, n = 30, design = "stratified",
                           coords = c("x", "y"), strata = "treat")
  expect_equal(nrow(p), 30L)
  expect_setequal(as.character(unique(p$treat)), levels(g$treat))
})

test_that("cell centres sit inside the trial", {
  g <- make_trial_grid(240, 108, cell_size = 9)
  expect_true(all(g$x > 0 & g$x < 240))
  expect_true(all(g$y > 0 & g$y < 108))
})

test_that("make_trial_grid validates its inputs", {
  expect_error(make_trial_grid(-1, 100), "positive")
  expect_error(make_trial_grid(240, 108, cell_size = 500), "no larger")
  expect_error(make_trial_grid(240, 108, treatments = "A"), "at least two")
  expect_error(make_trial_grid(240, 108, strip_width = -5), "positive")
})
