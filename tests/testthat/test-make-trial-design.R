test_that("a strip layout has the geometry the arguments describe", {
  d <- make_trial_design(c("A", "B", "C"), n_rep = 4, plot_width = 18,
                         plot_length = 200, cell_size = 9, seed = 1)
  a <- attr(d, "design")

  expect_equal(a$n_plots, 12L)
  expect_equal(a$trial_width, 12 * 18)
  expect_equal(a$trial_length, 200)
  expect_equal(nrow(d), max(d$row) * max(d$col))
  # Every treatment gets the same area
  expect_equal(length(unique(table(d$treat))), 1L)
})

test_that("a stacked layout halves the width and doubles the length", {
  strip <- attr(make_trial_design(c("A", "B", "C"), n_rep = 4,
                                  plot_width = 18, plot_length = 200,
                                  seed = 1), "design")
  stack <- attr(make_trial_design(c("A", "B", "C"), n_rep = 4,
                                  layout = "stack", plot_width = 18,
                                  plot_length = 200, gap = 10,
                                  seed = 1), "design")

  expect_equal(stack$trial_width, strip$trial_width / 2)
  expect_equal(stack$trial_length, 2 * 200 + 10)
  expect_equal(stack$n_plots, strip$n_plots)
})

test_that("the buffer between tiers carries no treatment", {
  s <- make_trial_design(c("A", "B", "C"), n_rep = 4, layout = "stack",
                         plot_width = 18, plot_length = 200, gap = 18,
                         cell_size = 9, seed = 1)
  expect_true(any(is.na(s$treat)))
  # Buffer cells have no plot or replicate either
  buf <- s[is.na(s$treat), ]
  expect_true(all(is.na(buf$plot)))
  expect_true(all(is.na(buf$rep)))
  # and they sit between the two tiers
  expect_true(all(is.na(s$tier[is.na(s$treat)])))
})

test_that("each replicate block contains every treatment exactly once", {
  d <- make_trial_design(c("A", "B", "C", "D"), n_rep = 3, seed = 2)
  book <- unique(d[!is.na(d$plot), c("rep", "plot", "treat")])
  for (r in unique(book$rep)) {
    expect_setequal(as.character(book$treat[book$rep == r]),
                    c("A", "B", "C", "D"))
  }
})

test_that("randomise = FALSE repeats one order, TRUE varies it", {
  syst <- make_trial_design(c("A", "B", "C"), n_rep = 4, randomise = FALSE)
  order_of <- function(d) {
    b <- unique(d[!is.na(d$plot), c("rep", "plot", "treat")])
    b <- b[order(b$plot), ]
    split(as.character(b$treat), b$rep)
  }
  os <- order_of(syst)
  expect_true(all(vapply(os, identical, logical(1), os[[1]])))

  # Over several seeds a randomised layout must produce more than one order
  seen <- unique(unlist(lapply(1:8, function(s) {
    paste(unlist(order_of(make_trial_design(c("A", "B", "C"), n_rep = 4,
                                            seed = s))), collapse = "")
  })))
  expect_gt(length(seen), 1L)
})

test_that("a design is reproducible from its seed", {
  a <- make_trial_design(c("A", "B", "C"), n_rep = 4, seed = 99)
  b <- make_trial_design(c("A", "B", "C"), n_rep = 4, seed = 99)
  expect_equal(a$treat, b$treat)
})

test_that("the design feeds place_point_samples directly", {
  d <- make_trial_design(c("A", "B", "C"), n_rep = 4, seed = 1)
  set.seed(1)
  p <- place_point_samples(d[!is.na(d$treat), ], n = 24, design = "stratified",
                           coords = c("x", "y"), strata = "treat")
  expect_equal(nrow(p), 24L)
  expect_setequal(as.character(unique(p$treat)), c("A", "B", "C"))
})

test_that("make_trial_design validates its inputs", {
  expect_error(make_trial_design("A"), "at least two")
  expect_error(make_trial_design(c("A", "B"), n_rep = 0), "at least 1")
  expect_error(make_trial_design(c("A", "B"), plot_width = -1), "positive")
  # A cell wider than a plot could never sit inside one treatment
  expect_error(make_trial_design(c("A", "B"), plot_width = 6, cell_size = 9),
               "wider than a plot")
  # Stacking needs an even number of plots to split into two tiers
  expect_error(make_trial_design(c("A", "B", "C"), n_rep = 3,
                                 layout = "stack"), "must be even")
})
