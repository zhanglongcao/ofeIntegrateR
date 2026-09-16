# A paddock with a ridge running across it and a sandy corner: two features in
# two covariates, so the structure is known and the zones can be checked.
ridge_paddock <- function(seed = 1, n_row = 24, n_col = 18) {
  set.seed(seed)
  g <- expand.grid(row = seq_len(n_row), col = seq_len(n_col))
  g$x <- g$col * 10
  g$y <- g$row * 10
  g$elevation <- 100 + 6 * exp(-((g$y - 120)^2) / 2000) +
    stats::rnorm(nrow(g), 0, 0.2)
  g$clay <- ifelse(g$x > 120 & g$y > 140, 18, 32) + stats::rnorm(nrow(g), 0, 1)
  g
}

two_block <- function(seed = 3, n = 20) {
  set.seed(seed)
  g <- expand.grid(row = seq_len(n), col = seq_len(n))
  g$x <- g$col * 10
  g$y <- g$row * 10
  g$ec <- ifelse(g$col > 10, 60, 40) + stats::rnorm(nrow(g), 0, 2)
  g
}

test_that("every skater zone is a single connected piece", {
  # The contract this method exists for. Checked on three different fields,
  # and at several k, because contiguity must not depend on the data.
  for (s in 1:3) {
    g <- ridge_paddock(seed = s)
    for (kk in c(2L, 4L, 6L)) {
      z <- partition_paddock(g, covariates = c("elevation", "clay"), k = kk)
      p <- attr(z, "partition")
      expect_equal(p$k, kk)
      expect_true(all(p$zones$patches == 1L),
                  info = paste("seed", s, "k", kk))
    }
  }
})

test_that("k-means is offered, and is shown to fragment", {
  # The comparison the default exists to win. If this ever stops fragmenting
  # the docs claiming skater is needed should be revisited, not the test.
  set.seed(11)
  g <- ridge_paddock()
  km <- partition_paddock(g, covariates = c("elevation", "clay"),
                          method = "kmeans", k = 4)
  sk <- partition_paddock(g, covariates = c("elevation", "clay"), k = 4)
  expect_gt(max(attr(km, "partition")$zones$patches), 1L)
  expect_true(all(attr(sk, "partition")$zones$patches == 1L))
  expect_equal(attr(km, "partition")$method, "kmeans")
})

test_that("a uniform paddock is reported as one zone", {
  set.seed(2)
  u <- expand.grid(row = 1:20, col = 1:20)
  u$x <- u$col * 10
  u$y <- u$row * 10
  u$elevation <- stats::rnorm(400, 100, 1)
  u$clay <- stats::rnorm(400, 30, 2)
  z <- partition_paddock(u, covariates = c("elevation", "clay"))
  expect_equal(attr(z, "partition")$k, 1L)
  expect_equal(nlevels(z$zone), 1L)
  expect_lt(attr(z, "partition")$table$r2[2], 0.10)
})

test_that("a planted two-block paddock is split in the right place", {
  g <- two_block()
  z <- partition_paddock(g, covariates = "ec")
  p <- attr(z, "partition")
  expect_equal(p$k, 2L)
  expect_true(all(p$zones$patches == 1L))
  # Zone 1 is the low-EC half by construction of the zone ordering.
  expect_lt(p$zones$ec[1], p$zones$ec[2])
  agree <- mean((as.integer(z$zone) == 2L) == (g$col > 10))
  expect_gt(agree, 0.97)
})

test_that("the criterion table is complete and drives the choice", {
  g <- ridge_paddock()
  z <- partition_paddock(g, covariates = c("elevation", "clay"))
  tab <- attr(z, "partition")$table
  expect_named(tab, c("k", "ssd", "r2", "gain", "ch"))
  expect_equal(tab$k, seq_len(nrow(tab)))
  expect_true(all(diff(tab$ssd) < 0))          # each split explains more
  expect_equal(tab$r2[1], 0)
  expect_equal(tab$gain[-1], diff(tab$r2))
  # The chosen k is the first whose successor is not worth having.
  k <- attr(z, "partition")$k
  expect_true(all(tab$gain[2:k] >= 0.05))
  if (k < nrow(tab)) expect_lt(tab$gain[k + 1], 0.05)

  # A stricter threshold can only give fewer zones.
  strict <- partition_paddock(g, covariates = c("elevation", "clay"),
                              min_gain = 0.20)
  expect_lte(attr(strict, "partition")$k, k)

  ch <- partition_paddock(g, covariates = c("elevation", "clay"),
                          criterion = "ch")
  expect_gte(attr(ch, "partition")$k, 1L)
})

test_that("min_cells is respected and limits what k can be asked for", {
  g <- ridge_paddock()
  z <- partition_paddock(g, covariates = c("elevation", "clay"), k = 5,
                         min_cells = 40)
  expect_true(all(attr(z, "partition")$zones$n >= 40))
  expect_equal(attr(z, "partition")$min_cells, 40L)
  expect_error(partition_paddock(g, covariates = "elevation", k = 20,
                                 min_cells = 50), "only 432 are available")
})

test_that("zones are numbered by the first covariate, so runs are comparable", {
  g <- ridge_paddock()
  a <- partition_paddock(g, covariates = c("elevation", "clay"), k = 3)
  b <- partition_paddock(g, covariates = c("elevation", "clay"), k = 3)
  expect_equal(as.integer(a$zone), as.integer(b$zone))
  expect_true(all(diff(attr(a, "partition")$zones$elevation) > 0))
})

test_that("cells missing a covariate are dropped, not guessed", {
  g <- ridge_paddock()
  g$clay[c(5, 40, 300)] <- NA
  z <- partition_paddock(g, covariates = c("elevation", "clay"), k = 3)
  expect_equal(nrow(z), nrow(g))
  expect_true(all(is.na(z$zone[c(5, 40, 300)])))
  expect_equal(attr(z, "partition")$n_dropped, 3L)
  expect_equal(attr(z, "partition")$n_cells, nrow(g) - 3L)
  expect_equal(sum(attr(z, "partition")$zones$n), nrow(g) - 3L)
})

test_that("irregular point layers are zoned through nearest neighbours", {
  set.seed(5)
  ir <- data.frame(x = stats::runif(300, 0, 200), y = stats::runif(300, 0, 300))
  ir$elevation <- 100 + 5 * (ir$y > 150) + stats::rnorm(300, 0, 0.5)
  z <- partition_paddock(ir, covariates = "elevation", x = "x", y = "y",
                         neighbours = "knn", n_neighbours = 6)
  p <- attr(z, "partition")
  expect_equal(p$k, 2L)
  expect_equal(p$neighbours, "knn")
  expect_true(all(p$zones$patches == 1L))
  expect_gt(mean((as.integer(z$zone) == 2L) == (ir$y > 150)), 0.95)
  # Grid neighbourhoods are not available without row/col.
  expect_error(partition_paddock(ir, covariates = "elevation", x = "x",
                                 y = "y", neighbours = "rook"),
               "needs `row` and `col`")
})

test_that("treatment coverage is reported and a gap is warned about", {
  g <- ridge_paddock()
  g$treat <- factor(c("A", "B", "C")[cut(g$col, 3, labels = FALSE)])
  expect_warning(z <- partition_paddock(g, covariates = c("elevation", "clay"),
                                        treat = "treat"),
                 "not estimable")
  summ <- attr(z, "partition")$zones
  expect_true("treatments" %in% names(summ))
  expect_true(any(summ$treatments < 3L))

  # Zones that do cover every treatment pass without complaint.
  suppressWarnings(z2 <- partition_paddock(g, covariates = "elevation", k = 2,
                                           treat = "treat"))
  expect_equal(attr(z2, "partition")$zones$treatments, c(3L, 3L))
})

test_that("input that cannot be zoned is refused with a reason", {
  g <- ridge_paddock()
  expect_error(partition_paddock(g, covariates = character(0)),
               "at least one column")
  expect_error(partition_paddock(g, covariates = "nope"), "not found in `data`")
  g$soil_type <- factor("sand")
  expect_error(partition_paddock(g, covariates = "soil_type"), "must be numeric")
  flat <- data.frame(row = 1:20, col = 1, x = 1:20, y = 1, z1 = 1, z2 = 2)
  expect_error(partition_paddock(flat, covariates = c("z1", "z2")),
               "nothing to partition")
  bare <- data.frame(elevation = stats::rnorm(20))
  expect_error(partition_paddock(bare, covariates = "elevation"),
               "Need coordinates")
})

test_that("a constant covariate is dropped rather than silently weighted", {
  g <- ridge_paddock()
  g$useless <- 7
  expect_warning(z <- partition_paddock(g, covariates = c("elevation",
                                                          "useless"), k = 2),
                 "Constant covariate")
  expect_equal(attr(z, "partition")$clustered_on, "elevation")
})

test_that("n_pc clusters on components instead of the covariates", {
  g <- ridge_paddock()
  g$elev2 <- g$elevation * 2 + stats::rnorm(nrow(g), 0, 0.05)
  z <- partition_paddock(g, covariates = c("elevation", "elev2", "clay"),
                         n_pc = 2, k = 3)
  expect_equal(attr(z, "partition")$n_pc, 2L)
  expect_true(all(attr(z, "partition")$zones$patches == 1L))
  # The covariate means are still reported on their original scale.
  expect_true(all(c("elevation", "elev2", "clay") %in%
                    names(attr(z, "partition")$zones)))
})

test_that("the zones feed the analysis functions", {
  g <- ridge_paddock()
  g$treat <- factor(rep(c("A", "B", "C"), length.out = nrow(g)))
  set.seed(9)
  g$yield <- 3 + 0.3 * as.integer(g$treat) + 0.2 * (g$elevation - 100) +
    stats::rnorm(nrow(g), 0, 0.3)
  suppressWarnings(z <- partition_paddock(g, covariates = c("elevation",
                                                            "clay"), k = 2))
  r <- adaptive_residual(z)
  expect_match(deparse(r), "dsum")
  fit <- fit_ofe(yield ~ zone + treat, data = z, residual = r)
  expect_s3_class(fit, "ofe_fit")
  expect_equal(length(fit$sections), 2L)
})

test_that("the spanning tree spans, and disconnected graphs are joined", {
  g <- two_block(n = 10)
  xy <- cbind(g$x, g$y)
  e <- ofeIntegrateR:::.cell_neighbours(xy, g$row, g$col, "rook")
  expect_equal(length(unique(ofeIntegrateR:::.components(e, nrow(g)))), 1L)
  w <- abs(g$ec[e[, 1]] - g$ec[e[, 2]])
  tr <- ofeIntegrateR:::.mst(e, w, nrow(g))
  expect_equal(nrow(tr), nrow(g) - 1L)
  expect_equal(length(unique(ofeIntegrateR:::.components(tr, nrow(g)))), 1L)

  # Two islands with no rook link between them are bridged, not rejected.
  split <- rbind(expand.grid(row = 1:4, col = 1:4),
                 expand.grid(row = 1:4, col = 7:10))
  sxy <- cbind(split$col * 10, split$row * 10)
  se <- ofeIntegrateR:::.cell_neighbours(sxy, split$row, split$col, "rook")
  expect_equal(length(unique(ofeIntegrateR:::.components(se, nrow(split)))), 2L)
  joined <- ofeIntegrateR:::.connect_components(se, sxy)
  expect_equal(length(unique(ofeIntegrateR:::.components(joined,
                                                         nrow(split)))), 1L)
})
