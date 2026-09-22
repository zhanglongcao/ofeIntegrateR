test_that("the table reports the trade-off a cell size actually makes", {
  sim <- simulate_yield_monitor(seed = 11)
  tab <- choose_cell_size(sim$cloud, response = "yield", treat = "treat")
  expect_named(tab, c("cell_size", "cells", "occupied", "coverage",
                      "median_obs", "pure", "usable"))
  expect_gt(nrow(tab), 3)
  expect_true(all(diff(tab$cell_size) > 0))

  # The trade-off itself: bigger cells hold more and cover more, but straddle
  # treatment boundaries more often. If these ever stop holding, the function
  # is measuring something other than what it claims.
  expect_true(all(diff(tab$median_obs) >= 0))
  expect_true(all(diff(tab$cells) < 0))
  expect_lt(tab$pure[nrow(tab)], tab$pure[1])
  expect_true(all(tab$occupied <= tab$cells))
  expect_true(all(tab$usable <= tab$occupied))
})

test_that("the recommendation follows the stated rule, and says what it is", {
  sim <- simulate_yield_monitor(seed = 11)
  tab <- choose_cell_size(sim$cloud, response = "yield", treat = "treat",
                          min_obs = 3, min_purity = 0.7)
  rec <- attr(tab, "recommended")
  ok <- tab$median_obs >= 3 & tab$pure >= 0.7
  expect_equal(rec, tab$cell_size[which(ok)[1]])
  expect_match(attr(tab, "rule"), "smallest size")

  # A stricter support floor moves the recommendation up.
  strict <- choose_cell_size(sim$cloud, response = "yield", treat = "treat",
                             min_obs = 5, min_purity = 0.7)
  expect_gt(attr(strict, "recommended"), rec)

  # Stricter still and there is no answer at all, because on this geometry --
  # 18 m treatment strips cut by 9 m passes -- a cell large enough to hold
  # eight observations already straddles a treatment boundary. The function
  # says so rather than returning the least bad size as though it qualified.
  expect_warning(none <- choose_cell_size(sim$cloud, response = "yield",
                                          treat = "treat", min_obs = 8,
                                          min_purity = 0.7),
                 "No size tried met both floors")
  expect_true(is.na(attr(none, "recommended")))
})

test_that("the recommendation is usable by the function it is for", {
  sim <- simulate_yield_monitor(seed = 11)
  tab <- choose_cell_size(sim$cloud, response = "yield", treat = "treat")
  g <- grid_dense_layer(sim$cloud, response = "yield", treat = "treat",
                        cell_size = attr(tab, "recommended"))
  expect_equal(nrow(g), tab$cells[tab$cell_size == attr(tab, "recommended")])
  expect_gte(stats::median(g$n_obs[g$n_obs > 0]), 3)
})

test_that("it works without a treatment column, and warns when nothing fits", {
  sim <- simulate_yield_monitor(seed = 11)
  tab <- choose_cell_size(sim$cloud, response = "yield")
  expect_true(all(is.na(tab$pure)))
  expect_false(is.na(attr(tab, "recommended")))

  expect_warning(
    choose_cell_size(sim$cloud, response = "yield", treat = "treat",
                     sizes = c(30, 40), min_purity = 0.95),
    "No size tried met both floors")
})

test_that("bad input is refused", {
  sim <- simulate_yield_monitor(seed = 11)
  expect_error(choose_cell_size(sim$cloud, response = "nope"),
               "not found in `data`")
  expect_error(choose_cell_size(sim$cloud, response = "yield", sizes = -1),
               "must contain a positive value")
})
