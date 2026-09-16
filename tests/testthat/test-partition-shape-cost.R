# The help page and README quote how much rectangles cost against skater on
# differently shaped features. Those numbers are a claim about behaviour, so
# they are checked here rather than left to rot.
shaped <- function(kind, seed = 1) {
  set.seed(seed)
  g <- expand.grid(row = 1:24, col = 1:24)
  g$x <- g$col * 10
  g$y <- g$row * 10
  g$v <- switch(
    kind,
    blocks   = ifelse(g$x > 120, 10, 0) + ifelse(g$y > 120, 5, 0),
    diagonal = 10 * ((g$x + g$y) > 240),
    contour  = 10 * (sqrt((g$x - 120)^2 + (g$y - 120)^2) < 80),
    gradient = 0.04 * g$y)
  g$v <- g$v + stats::rnorm(nrow(g), 0, 0.8)
  g
}
r2_at <- function(g, method, k) {
  attr(partition_paddock(g, covariates = "v", k = k, method = method),
       "partition")$table$r2[k]
}

test_that("rectangles cost nothing on axis-aligned or gradient structure", {
  for (kind in c("blocks", "gradient")) {
    rect <- r2_at(shaped(kind), "rectangle", 2L)
    skat <- r2_at(shaped(kind), "skater", 2L)
    expect_gt(rect, 0.6)
    expect_gt(rect, skat - 0.02)
  }
})

test_that("rectangles cost heavily on a diagonal or curved boundary", {
  for (kind in c("diagonal", "contour")) {
    rect <- r2_at(shaped(kind), "rectangle", 2L)
    skat <- r2_at(shaped(kind), "skater", 2L)
    expect_gt(skat, 0.9)
    expect_lt(rect, 0.4)
    # ... and raising k lets rectangles approximate it in steps, as documented.
    expect_gt(r2_at(shaped(kind), "rectangle", 6L), rect)
  }
})
