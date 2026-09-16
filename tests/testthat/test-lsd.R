balanced <- function(k = 6, r = 8, seed = 20,
                     mu = c(10, 10.6, 11.2, 12.5, 12.7, 15)) {
  set.seed(seed)
  d <- data.frame(treat = factor(rep(LETTERS[seq_len(k)], each = r)),
                  row = rep(seq_len(r), k), col = rep(seq_len(k), each = r))
  d$y <- mu[as.integer(d$treat)] + stats::rnorm(k * r, 0, 1)
  d
}

test_that("compact_letters reproduces the textbook displays", {
  three <- function(p) data.frame(level1 = c("A", "A", "B"),
                                  level2 = c("B", "C", "C"), p.value = p)
  mu <- c(A = 1, B = 2, C = 3)

  # A and B alike, both apart from C
  expect_equal(compact_letters(three(c(0.40, 0.001, 0.002)), means = mu),
               c(C = "a", B = "b", A = "b"))
  # A chain: the middle mean shares a letter with each end, the ends do not
  expect_equal(compact_letters(three(c(0.20, 0.001, 0.30)), means = mu),
               c(C = "a", B = "ab", A = "b"))
  expect_equal(compact_letters(three(c(0.001, 0.001, 0.001)), means = mu),
               c(C = "a", B = "b", A = "c"))
  expect_equal(compact_letters(three(c(0.9, 0.9, 0.9)), means = mu),
               c(C = "a", B = "a", A = "a"))
  # Smaller is better
  expect_equal(compact_letters(three(c(0.40, 0.001, 0.002)), means = mu,
                               decreasing = FALSE),
               c(A = "a", B = "a", C = "b"))
})

test_that("compact_letters takes a matrix, a contrast column, or level columns", {
  pm <- matrix(c(NA, 0.4, 0.001, 0.4, NA, 0.002, 0.001, 0.002, NA), 3, 3,
               dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_equal(compact_letters(pm, means = c(A = 1, B = 2, C = 3)),
               c(C = "a", B = "b", A = "b"))

  by_contrast <- data.frame(contrast = c("B - A", "C - A", "C - B"),
                            p.value = c(0.4, 0.001, 0.002))
  expect_equal(compact_letters(by_contrast, means = c(A = 1, B = 2, C = 3)),
               c(C = "a", B = "b", A = "b"))

  # Without means the levels keep the order they appear in.
  expect_named(compact_letters(by_contrast), c("A", "B", "C"))
})

test_that("compact_letters refuses input it cannot letter honestly", {
  cmp <- data.frame(level1 = "A", level2 = "B", p.value = 0.3)
  expect_error(compact_letters(cmp, p_col = "prob"), "no column `prob`")
  expect_error(compact_letters(cmp, means = c(1, 2)), "must be named")
  expect_error(compact_letters(cmp, means = c(A = 1)), "missing level")
  expect_error(compact_letters(matrix(0.5, 2, 2)), "needs dimnames")
  # A pair that was never compared must not pass silently as "different".
  gap <- data.frame(level1 = c("A", "A"), level2 = c("B", "C"),
                    p.value = c(0.001, 0.001))
  expect_warning(compact_letters(gap), "no p-value")
})

test_that("levels share a letter exactly when they are not separable", {
  # The defining property of a compact letter display, checked against the
  # comparisons that produced it rather than against a stored answer.
  fit <- fit_ofe(y ~ treat, data = balanced(), residual = ~ id(row):id(col))
  tab <- ofe_lsd(fit, "treat")
  cmp <- attr(tab, "comparisons")
  grp <- stats::setNames(tab$group, tab$treat)
  for (i in seq_len(nrow(cmp))) {
    shared <- length(intersect(
      strsplit(grp[[cmp$level1[i]]], "")[[1]],
      strsplit(grp[[cmp$level2[i]]], "")[[1]])) > 0
    expect_equal(shared, cmp$p.value[i] > 0.05,
                 info = cmp$contrast[i])
  }
})

test_that("the reported LSD is t times the average SED", {
  fit <- fit_ofe(y ~ treat, data = balanced(), residual = ~ id(row):id(col))
  tab <- ofe_lsd(fit, "treat")
  l <- attr(tab, "lsd")
  cmp <- attr(tab, "comparisons")
  expect_equal(l$average_sed, mean(cmp$std.error))
  expect_equal(l$lsd, stats::qt(0.975, l$df) * l$average_sed)
  expect_equal(l$df, fit$n - fit$p)

  # A tighter alpha gives a wider LSD and never more letters.
  t01 <- ofe_lsd(fit, "treat", alpha = 0.01)
  expect_gt(attr(t01, "lsd")$lsd, l$lsd)
  expect_lte(length(unique(unlist(strsplit(t01$group, "")))),
             length(unique(unlist(strsplit(tab$group, "")))))
})

test_that("ofe_lsd output is ordered, labelled and complete", {
  fit <- fit_ofe(y ~ treat, data = balanced(), residual = ~ id(row):id(col))
  tab <- ofe_lsd(fit, "treat")
  expect_named(tab, c("treat", "estimate", "std.error", "lower", "upper",
                      "group"))
  expect_equal(nrow(tab), 6L)
  expect_true(all(diff(tab$estimate) < 0))            # sorted, best first
  expect_equal(tab$group[1], "a")
  expect_true(all(tab$lower < tab$estimate & tab$estimate < tab$upper))

  unsorted <- ofe_lsd(fit, "treat", sort = FALSE)
  expect_equal(unsorted$treat, LETTERS[1:6])

  # The means must be the ones ofe_means() reports.
  mu <- ofe_means(fit, "treat")
  expect_equal(tab$estimate[match(mu$treat, tab$treat)], mu$estimate)
})

test_that("a multiplicity adjustment can only merge letters, never split them", {
  fit <- fit_ofe(y ~ treat, data = balanced(), residual = ~ id(row):id(col))
  none <- ofe_lsd(fit, "treat")
  for (adj in c("tukey", "bonferroni", "holm", "sidak", "BH")) {
    adjusted <- ofe_lsd(fit, "treat", adjust = adj)
    expect_gte(min(attr(adjusted, "comparisons")$p.value),
               min(attr(none, "comparisons")$p.value))
    expect_lte(length(unique(unlist(strsplit(adjusted$group, "")))),
               length(unique(unlist(strsplit(none$group, "")))))
  }
  expect_error(ofe_lsd(fit, "treat", adjust = "tukey", use = "lsd"),
               "single average LSD")
  expect_error(ofe_lsd(fit, "treat", adjust = "nonsense"), "should be one of")
})

test_that("use = 'lsd' judges every pair against the same difference", {
  fit <- fit_ofe(y ~ treat, data = balanced(), residual = ~ id(row):id(col))
  tab <- ofe_lsd(fit, "treat", use = "lsd")
  l <- attr(tab, "lsd")
  cmp <- attr(tab, "comparisons")
  grp <- stats::setNames(tab$group, tab$treat)
  for (i in seq_len(nrow(cmp))) {
    shared <- length(intersect(
      strsplit(grp[[cmp$level1[i]]], "")[[1]],
      strsplit(grp[[cmp$level2[i]]], "")[[1]])) > 0
    expect_equal(shared, abs(cmp$estimate[i]) <= l$lsd, info = cmp$contrast[i])
  }
})

test_that("by = gives one lettering per level, and never letters across them", {
  set.seed(4)
  d <- balanced(k = 3, r = 20, mu = c(10, 11, 12))
  d$zone <- factor(ifelse(d$row <= 10, "north", "south"))
  # The treatment effect is doubled in the south, so the letterings differ.
  d$y <- d$y + ifelse(d$zone == "south", as.numeric(d$treat) * 2, 0)
  fit <- fit_ofe(y ~ zone + zone:treat, data = d, residual = ~ id(row):id(col))

  tab <- ofe_lsd(fit, "treat", by = "zone")
  expect_equal(names(tab)[1:2], c("zone", "treat"))
  expect_equal(nrow(tab), 6L)
  expect_equal(sort(unique(tab$zone)), c("north", "south"))
  expect_true(all(tapply(tab$group, tab$zone, function(g) g[1]) == "a"))

  # Comparisons stay inside a zone: three pairs each, none crossing.
  cmp <- attr(tab, "comparisons")
  expect_equal(nrow(cmp), 6L)
  expect_equal(as.vector(table(cmp$zone)), c(3L, 3L))
  # One LSD per zone, not one overall.
  expect_equal(nrow(attr(tab, "lsd")), 2L)

  # The by-means match the same model's means for that zone.
  mu <- ofe_means(fit, "treat", by = "zone")
  key <- paste(tab$zone, tab$treat)
  expect_equal(tab$estimate[match(paste(mu$zone, mu$treat), key)], mu$estimate)
})

test_that("ofe_means gains by and adjust without changing what it meant", {
  fit <- fit_ofe(y ~ treat, data = balanced(), residual = ~ id(row):id(col))
  pw <- ofe_means(fit, "treat", pairwise = TRUE)
  expect_true(all(c("level1", "level2", "contrast") %in% names(pw)))
  expect_equal(nrow(pw), 15L)

  tk <- ofe_means(fit, "treat", pairwise = TRUE, adjust = "tukey")
  expect_true(all(tk$p.value >= pw$p.value))
  expect_equal(tk$estimate, pw$estimate)

  expect_error(ofe_means(fit, "nope"), "not a column")
  expect_error(ofe_lsd(fit, "treat", by = "nope"), "not a column")
})

test_that("ofe_lsd agrees with agricolae on a balanced design", {
  skip_if_not_installed("agricolae")
  d <- balanced()
  ref <- agricolae::LSD.test(stats::aov(y ~ treat, data = d), "treat",
                             console = FALSE)
  fit <- fit_ofe(y ~ treat, data = d, residual = ~ id(row):id(col))
  mine <- ofe_lsd(fit, "treat", use = "lsd")
  expect_equal(attr(mine, "lsd")$lsd, ref$statistics$LSD, tolerance = 1e-8)
  expect_equal(mine$estimate,
               ref$groups$y[match(mine$treat, rownames(ref$groups))],
               tolerance = 1e-8)
  expect_equal(mine$group,
               as.character(ref$groups$groups)[match(mine$treat,
                                                     rownames(ref$groups))])

  hsd <- agricolae::HSD.test(stats::aov(y ~ treat, data = d), "treat",
                             console = FALSE)
  tk <- ofe_lsd(fit, "treat", adjust = "tukey")
  expect_equal(tk$group,
               as.character(hsd$groups$groups)[match(tk$treat,
                                                     rownames(hsd$groups))])
})
