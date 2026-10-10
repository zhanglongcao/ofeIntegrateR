test_that("place_count_stations lays out stations inside the buffers", {
  d <- make_trial_design(c("A", "B"), n_rep = 3, plot_width = 12,
                         plot_length = 250, seed = 1)
  st <- place_count_stations(d, n_stations = 8, n_quadrats = 2,
                             end_buffer = 20, side_buffer = 2, seed = 1)
  expect_s3_class(st, "ofe_stations")
  expect_equal(nrow(st), 6 * 8 * 2)
  expect_true(all(st$along >= 20 & st$along <= 230))
  pl <- attr(st, "plots")
  expect_equal(pl$xmax - pl$xmin, rep(12, 6))
  for (k in seq_len(nrow(pl))) {
    xs <- st$x[st$plot == pl$plot[k]]
    expect_true(all(xs >= pl$xmin[k] + 2 & xs <= pl$xmax[k] - 2))
  }
})

test_that("systematic stations sit at the same positions in every plot", {
  d <- make_trial_design(c("A", "B"), n_rep = 2, plot_width = 12,
                         plot_length = 200, seed = 1)
  st <- place_count_stations(d, n_stations = 6, n_quadrats = 1, seed = 2)
  pos <- split(st$along, st$plot)
  for (k in 2:length(pos)) expect_equal(pos[[k]], pos[[1]], tolerance = 1e-8)
  expect_equal(diff(pos[[1]]), rep((200 - 40) / 6, 5), tolerance = 1e-8)
})

test_that("stacked layouts and bad buffers are handled", {
  d <- make_trial_design(c("A", "B"), n_rep = 2, layout = "stack",
                         plot_width = 12, plot_length = 150, gap = 10, seed = 1)
  st <- place_count_stations(d, 5, 2, seed = 1)
  pl <- attr(st, "plots")
  expect_equal(sort(unique(pl$ymin)), c(0, 160))
  expect_error(place_count_stations(d, 5, 2, end_buffer = 80), "no room")
  expect_error(place_count_stations(d, 5, 2, side_buffer = 6), "no room")
  expect_error(place_count_stations(d, 0, 2), "positive")
})

test_that("count_precision falls with more stations and is reproducible", {
  a <- count_precision(c(4, 16), n_quadrats = 2, quadrat_area = 0.25,
                       crop = "wheat", n_sim = 60, seed = 3)
  b <- count_precision(c(4, 16), n_quadrats = 2, quadrat_area = 0.25,
                       crop = "wheat", n_sim = 60, seed = 3)
  expect_equal(a, b)
  expect_true(a$margin95[2] < a$margin95[1])
  expect_lt(abs(mean(a$bias)), 0.05)
  expect_error(count_precision(5, crop = "custom", n_sim = 5), "supply")
})

test_that("count_sample_size returns the smallest adequate design", {
  r <- count_sample_size(0.5, max_stations = 10, crop = "wheat",
                         quadrat_area = 0.25, n_quadrats = 1, n_sim = 40, seed = 1)
  expect_equal(nrow(r), 1L)
  expect_true(r$margin95 <= 0.5)
  expect_warning(count_sample_size(0.001, max_stations = 4, crop = "weeds",
                                   n_sim = 20, seed = 1), "No design")
})

test_that("analyse_strip_counts equals the paired t-test on log strip means", {
  set.seed(4)
  s <- data.frame(strip = 1:10, rep = rep(1:5, each = 2),
                  treat = rep(c("Ctl", "Trt"), 5))
  s$dens <- 100 * ifelse(s$treat == "Trt", 1.15, 1) * exp(rnorm(10, 0, 0.1))
  q <- s[rep(1:10, each = 12), ]
  q$count <- rpois(nrow(q), q$dens * 0.5)
  fit <- analyse_strip_counts(q, "count", "strip", "treat", block = "rep",
                              area = 0.5)
  r <- fit$strips$response
  dd <- r[fit$strips$treat == "Trt"] - r[fit$strips$treat == "Ctl"]
  expect_equal(fit$ratios$log_ratio, mean(dd))
  expect_equal(fit$ratios$p_value, t.test(dd)$p.value)
  expect_equal(nrow(fit$strips), 10L)
  expect_output(print(fit), "Ratio to Ctl")
})

test_that("analyse_strip_counts uses the baseline and checks its input", {
  set.seed(5)
  q <- data.frame(strip = rep(1:8, each = 6), rep = rep(rep(1:4, each = 2), each = 6),
                  treat = factor(rep(rep(c("Unsprayed", "Sprayed"), 4), each = 6),
                                 levels = c("Unsprayed", "Sprayed")))
  q$before <- rpois(nrow(q), 8)
  q$after <- rpois(nrow(q), q$before * ifelse(q$treat == "Sprayed", 0.5, 1))
  fit <- analyse_strip_counts(q, "after", "strip", "treat", "rep", area = 1,
                              baseline = "before")
  expect_true(fit$baseline)
  expect_equal(fit$ratios$reference, "Unsprayed")
  expect_true(fit$ratios$ratio < 1)
  expect_error(analyse_strip_counts(q, "nope", "strip", "treat"), "not found")
  q2 <- q; q2$treat[1] <- "Sprayed"   # strip 1 is Unsprayed
  expect_error(analyse_strip_counts(q2, "after", "strip", "treat"), "more than one")
})
