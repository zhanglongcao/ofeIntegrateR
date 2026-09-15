grid20 <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 5,
                             seed = 1)$grid

test_that("every design returns exactly n distinct locations", {
  set.seed(1)
  for (d in c("random", "grid", "stratified", "nested")) {
    p <- place_point_samples(grid20, n = 12, design = d)
    expect_equal(nrow(p), 12L, info = d)
    expect_equal(anyDuplicated(paste(p$row, p$col)), 0L, info = d)
    expect_equal(attr(p, "design"), d)
  }
})

test_that("stratified placement represents every stratum", {
  set.seed(2)
  p <- place_point_samples(grid20, n = 12, design = "stratified",
                           strata = "treat")
  expect_setequal(as.character(unique(p$treat)), levels(grid20$treat))
  # Balance should be much better than simple random placement
  expect_lt(diff(range(table(p$treat))), 4)
})

test_that("grid placement spreads more evenly than nested clustering", {
  set.seed(3)
  spread <- function(p) mean(dist(cbind(p$col, p$row)))
  g <- place_point_samples(grid20, n = 12, design = "grid")
  n <- place_point_samples(grid20, n = 12, design = "nested")
  expect_gt(spread(g), spread(n))
})

test_that("place_point_samples validates its inputs", {
  expect_error(place_point_samples(grid20, n = 12, coords = c("lon", "lat")),
               "coords")
  expect_error(place_point_samples(grid20, n = 12, design = "stratified",
                                   strata = "nope"), "not found")
  expect_error(place_point_samples(grid20, n = 0), "positive integer")
  expect_error(place_point_samples(grid20, n = nrow(grid20) + 1), "exceeds")
})
