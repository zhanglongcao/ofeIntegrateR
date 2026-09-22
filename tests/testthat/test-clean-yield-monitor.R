dirty <- function(seed = 1) simulate_yield_monitor(defects = TRUE, seed = seed)

test_that("the simulator labels the points it damages, and only those", {
  s <- dirty()
  d <- s$cloud
  expect_true(all(c("pass", "seq", "speed", "defect") %in% names(d)))
  expect_setequal(stats::na.omit(unique(d$defect)),
                  c("pass_start", "pass_end", "outlier"))
  # Pass starts read low and pass ends read high: that is the grain-flow lag,
  # and it is what makes these defects survive averaging.
  clean_med <- stats::median(d$yield[is.na(d$defect)])
  expect_lt(stats::median(d$yield[which(d$defect == "pass_start")]), clean_med)
  expect_gt(stats::median(d$yield[which(d$defect == "pass_end")]), clean_med)
  # Off by default, so nothing that existed before this changed.
  expect_true(all(is.na(simulate_yield_monitor(seed = 1)$cloud$defect)))
  expect_error(simulate_yield_monitor(defects = list(nope = 1)), "Unknown")
})

test_that("cleaning removes the damaged points and keeps the sound ones", {
  s <- dirty()
  d <- s$cloud
  f <- clean_yield_monitor(d, pass = "pass", order = "seq", speed = "speed",
                           pass_trim = 12, local_mad = 3, drop = FALSE)
  truth <- ifelse(is.na(d$defect), "clean", d$defect)
  recall <- function(k) mean(!f$.keep[truth == k])
  # The numbers quoted in the help page. They are a claim about behaviour, so
  # they are pinned here rather than left to drift.
  expect_gt(recall("pass_start"), 0.90)
  expect_gt(recall("pass_end"), 0.85)
  expect_lt(mean(!f$.keep[truth == "clean"]), 0.03)
})

test_that("a dropout is caught by min_yield, not by a dispersion rule", {
  s <- dirty()
  d <- s$cloud
  d$yield <- d$yield + 3                       # put it on a physical scale
  d$yield[which(s$cloud$defect == "outlier" & s$cloud$yield == 0)] <- 0
  zeros <- which(d$yield == 0)
  expect_gt(length(zeros), 5)

  without <- clean_yield_monitor(d, pass = "pass", speed = "speed",
                                 drop = FALSE)
  with_min <- clean_yield_monitor(d, pass = "pass", speed = "speed",
                                  min_yield = 0, drop = FALSE)
  # A zero sits well inside four MADs of the median, so the statistical rule
  # alone leaves it in. That is the whole reason the rule exists.
  expect_lt(mean(!without$.keep[zeros]), 0.5)
  expect_true(all(!with_min$.keep[zeros]))
  expect_true(all(with_min$.reason[zeros] == "zero"))
})

test_that("each removal is attributed to exactly one rule, and they tally", {
  s <- dirty()
  f <- clean_yield_monitor(s$cloud, pass = "pass", speed = "speed",
                           pass_trim = 10, drop = FALSE)
  tally <- attr(f, "cleaning")
  expect_equal(sum(tally$removed), nrow(s$cloud))
  expect_equal(tally$removed[tally$rule == "kept"], sum(f$.keep))
  for (r in setdiff(tally$rule, "kept")) {
    expect_equal(tally$removed[tally$rule == r], sum(f$.reason == r))
  }
  expect_equal(attr(f, "n_input"), nrow(s$cloud))
  expect_true(is.factor(f$.reason))
})

test_that("drop = TRUE returns the survivors and nothing else", {
  s <- dirty()
  kept <- clean_yield_monitor(s$cloud, pass = "pass", speed = "speed",
                              pass_trim = 10)
  audit <- clean_yield_monitor(s$cloud, pass = "pass", speed = "speed",
                               pass_trim = 10, drop = FALSE)
  expect_equal(nrow(kept), sum(audit$.keep))
  expect_false(".keep" %in% names(kept))
  expect_equal(names(kept), names(s$cloud))
  expect_equal(attr(kept, "cleaning"), attr(audit, "cleaning"))
})

test_that("passes are recovered from recording order and from geometry alike", {
  s <- dirty()
  d <- s$cloud
  from_order <- ofeIntegrateR:::.derive_passes(d$x, d$y, order = d$seq)
  from_geom <- ofeIntegrateR:::.derive_passes(d$x, d$y)
  # The obvious measure -- the share of points in each derived group's
  # majority true group -- scores a perfect 1 when every pass is shattered
  # into fragments, because it rewards purity and ignores completeness. The
  # count has to be checked too, or over-splitting passes unnoticed.
  expect_equal(length(unique(from_order)), length(unique(d$pass)))
  expect_equal(length(unique(from_geom)), length(unique(d$pass)))

  misassigned <- function(p) {
    tab <- table(d$pass, p)
    sum(tab) - sum(apply(tab, 1, max))
  }
  # Geometry recovers the passes exactly; recording order misplaces a handful
  # of points, and only ever at the joins between passes -- which is harmless,
  # because those are the points pass_trim removes whichever pass they land in.
  expect_equal(misassigned(from_geom), 0)
  expect_lt(misassigned(from_order) / nrow(d), 0.01)
  if (misassigned(from_order) > 0) {
    best <- apply(table(d$pass, from_order), 1, which.max)
    wrong <- which(from_order != best[as.character(d$pass)])
    from_end <- stats::ave(d$x, d$pass,
                           FUN = function(z) pmin(z - min(z), max(z) - z))
    expect_lt(max(from_end[wrong]), 12)
  }
})

test_that("trimming pass ends works without being told where the passes are", {
  s <- dirty()
  d <- s$cloud
  told <- clean_yield_monitor(d, pass = "pass", pass_trim = 12, drop = FALSE)
  guessed <- clean_yield_monitor(d, order = "seq", pass_trim = 12,
                                 drop = FALSE)
  expect_gt(mean(told$.keep == guessed$.keep), 0.95)
})

test_that("the edge rule removes a boundary band and nothing inside it", {
  set.seed(1)
  g <- expand.grid(x = seq(0, 100, 5), y = seq(0, 60, 5))
  g$yield <- 3 + stats::rnorm(nrow(g), 0, 0.1)
  f <- clean_yield_monitor(g, edge_buffer = 6, drop = FALSE)
  removed <- !f$.keep
  expect_true(all(f$.reason[removed] == "edge"))
  # Everything removed is near an edge; everything kept is not.
  near <- g$x <= 5 | g$x >= 95 | g$y <= 5 | g$y >= 55
  expect_true(all(near[removed]))
  expect_false(any(removed[g$x > 20 & g$x < 80 & g$y > 20 & g$y < 40]))
})

test_that("rules can be turned off and bad arguments are refused", {
  s <- dirty()
  none <- clean_yield_monitor(s$cloud, outlier_mad = Inf)
  expect_equal(nrow(none), nrow(s$cloud))

  expect_error(clean_yield_monitor(s$cloud, yield = "nope"),
               "not found in `data`")
  expect_error(clean_yield_monitor(s$cloud, yield_range = 3),
               "two increasing numbers")
  expect_error(clean_yield_monitor(s$cloud, yield_range = c(9, 1)),
               "two increasing numbers")
})

test_that("cleaning feeds the gridding step it exists to protect", {
  s <- dirty()
  kept <- clean_yield_monitor(s$cloud, pass = "pass", order = "seq",
                              speed = "speed", pass_trim = 12)
  g <- grid_dense_layer(kept, response = "yield", treat = "treat",
                        cell_size = 9)
  expect_true(all(c("row", "col", "yield", "n_obs") %in% names(g)))
  expect_gt(sum(!is.na(g$yield)), 50)
})
