make_grid <- function(n_row = 20, n_col = 12, seed = 1) {
  set.seed(seed)
  d <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))
  d$treat <- factor(c("A", "B", "C")[cut(d$col, 3, labels = FALSE)])
  d$rep <- factor(cut(d$row, 4, labels = FALSE))
  d
}

# Draw from a separable AR1 x AR1 field, which is the process fit_ofe claims
# to estimate; anything less would test the optimiser against itself.
ar1_field <- function(n_row, n_col, phi_r, phi_c, sd = 1) {
  Rr <- phi_r^abs(outer(seq_len(n_row), seq_len(n_row), "-"))
  Rc <- phi_c^abs(outer(seq_len(n_col), seq_len(n_col), "-"))
  L <- kronecker(chol(Rr), chol(Rc))
  sd * as.numeric(crossprod(L, stats::rnorm(n_row * n_col)))
}

test_that("an id():id() residual reproduces ordinary least squares exactly", {
  d <- make_grid()
  d$y <- as.numeric(d$treat) + stats::rnorm(nrow(d))
  fit <- fit_ofe(y ~ treat, data = d, residual = ~ id(row):id(col))
  ref <- stats::lm(y ~ treat, data = d)
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(vcov(fit)))),
               unname(sqrt(diag(vcov(ref)))), tolerance = 1e-6)
  expect_equal(fit$sigma2, summary(ref)$sigma^2, tolerance = 1e-8)
})

test_that("AR1 correlations and treatment effects are recovered", {
  set.seed(42)
  n_row <- 24; n_col <- 12
  d <- make_grid(n_row, n_col, seed = 42)
  d$y <- c(0, 1, 2)[as.integer(d$treat)] +
    ar1_field(n_row, n_col, phi_r = 0.7, phi_c = 0.4, sd = 1)
  fit <- fit_ofe(y ~ treat, data = d)
  expect_true(fit$converged)
  vc <- stats::setNames(fit$varcomp$estimate, fit$varcomp$component)
  expect_equal(unname(vc["R!row!cor"]), 0.7, tolerance = 0.25)
  expect_equal(unname(vc["R!col!cor"]), 0.4, tolerance = 0.3)
  est <- coef(fit)
  expect_equal(unname(est["treatB"]), 1, tolerance = 0.6)
  expect_equal(unname(est["treatC"]), 2, tolerance = 0.6)
})

test_that("a random term adds a variance component and shrinks nothing else", {
  d <- make_grid(24, 12, seed = 3)
  set.seed(3)
  d$y <- as.numeric(d$treat) + stats::rnorm(nlevels(d$rep), 0, 1)[d$rep] +
    stats::rnorm(nrow(d), 0, 0.5)
  fit <- fit_ofe(y ~ treat, data = d, random = ~ rep,
                 residual = ~ id(row):id(col))
  expect_true("rep" %in% fit$varcomp$component)
  expect_gt(fit$varcomp$variance[fit$varcomp$component == "rep"], 0)
  expect_error(fit_ofe(y ~ treat, data = d, random = ~ 1), "no terms")
})

test_that("dsum() gives each section its own parameters", {
  d <- make_grid(24, 12, seed = 5)
  d$zone <- factor(ifelse(d$row <= 12, "1", "2"))
  set.seed(5)
  d$y <- as.numeric(d$treat) + stats::rnorm(nrow(d), 0, ifelse(d$zone == "1", 0.5, 2))
  fit <- fit_ofe(y ~ treat, data = d,
                 residual = ~ dsum(~ id(row):id(col) | zone))
  # Section 2 was simulated with 4x the standard deviation, so 16x the variance.
  v <- fit$varcomp$estimate[fit$varcomp$component == "zone!2!var"]
  expect_gt(v, 8)
  expect_lt(v, 32)
})

test_that("missing responses are dropped and reported", {
  d <- make_grid()
  d$y <- as.numeric(d$treat) + stats::rnorm(nrow(d))
  d$y[c(3, 17, 44)] <- NA
  fit <- fit_ofe(y ~ treat, data = d, residual = ~ id(row):id(col))
  expect_equal(fit$n, nrow(d) - 3L)
  expect_equal(fit$n_dropped, 3L)
  expect_equal(nobs(fit), nrow(d) - 3L)
})

test_that("wald_tests reports one row per fixed term with the right df", {
  d <- make_grid()
  set.seed(9)
  d$x <- stats::rnorm(nrow(d))
  d$y <- as.numeric(d$treat) + d$x + stats::rnorm(nrow(d))
  fit <- fit_ofe(y ~ treat + x, data = d, residual = ~ id(row):id(col))
  w <- wald_tests(fit)
  expect_equal(w$term, c("treat", "x"))
  expect_equal(w$df, c(2L, 1L))
  expect_true(all(w$p.value >= 0 & w$p.value <= 1))
  expect_lt(w$p.value[w$term == "treat"], 0.01)
})

test_that("ofe_means agrees with the coefficients it is built from", {
  d <- make_grid()
  set.seed(11)
  d$y <- as.numeric(d$treat) + stats::rnorm(nrow(d))
  fit <- fit_ofe(y ~ treat, data = d, residual = ~ id(row):id(col))
  mu <- ofe_means(fit, "treat")
  expect_equal(names(mu)[1], "treat")
  expect_equal(nrow(mu), 3L)
  expect_equal(mu$estimate[2] - mu$estimate[1],
               unname(coef(fit)["treatB"]), tolerance = 1e-8)
  expect_true(all(mu$lower < mu$estimate & mu$estimate < mu$upper))

  pw <- ofe_means(fit, "treat", pairwise = TRUE)
  expect_equal(nrow(pw), 3L)
  expect_equal(pw$estimate[pw$contrast == "B - A"],
               unname(coef(fit)["treatB"]), tolerance = 1e-8)
  expect_equal(pw$std.error[pw$contrast == "B - A"],
               unname(sqrt(diag(vcov(fit)))["treatB"]), tolerance = 1e-8)
})

test_that("extract_fixed_effects handles an ofe_fit", {
  d <- make_grid()
  set.seed(13)
  d$y <- as.numeric(d$treat) + stats::rnorm(nrow(d))
  fit <- fit_ofe(y ~ treat, data = d, residual = ~ id(row):id(col))
  fe <- extract_fixed_effects(fit)
  expect_named(fe, c("term", "estimate", "se"))
  expect_equal(fe$estimate, unname(coef(fit)))
})

test_that("bad input is refused with a useful message", {
  d <- make_grid()
  d$y <- as.numeric(d$treat) + stats::rnorm(nrow(d))
  expect_error(fit_ofe(y ~ nope, data = d), "not found in `data`")
  expect_error(fit_ofe(~ treat, data = d), "two-sided formula")
  expect_error(fit_ofe(y ~ treat, data = d,
                       control = ofe_control(max_n = 10)), "max_n")
  d$flat <- factor("only")
  expect_error(fit_ofe(y ~ flat, data = d), "single level")
})

test_that("summary and print methods run and carry the fit", {
  d <- make_grid()
  set.seed(17)
  d$y <- as.numeric(d$treat) + stats::rnorm(nrow(d))
  fit <- fit_ofe(y ~ treat, data = d, residual = ~ id(row):id(col))
  s <- summary(fit)
  expect_s3_class(s, "summary.ofe_fit")
  expect_equal(nrow(s$coefficients), 3L)
  expect_output(print(fit), "Spatial mixed model")
  expect_output(print(s), "Fixed effects")
  expect_equal(length(fitted(fit)), fit$n)
  expect_equal(length(residuals(fit)), fit$n)
  expect_s3_class(logLik(fit), "logLik")
})

test_that("fit_integrated_kriged defaults to the licence-free engine", {
  sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 25, seed = 1)
  kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
  fit <- fit_integrated_kriged(kr, response = "dense_response", treat = "treat",
                               covariate = "point_obs_kriged")
  expect_s3_class(fit, "ofe_fit")
  expect_true("point_obs_kriged" %in% names(coef(fit)))
})
