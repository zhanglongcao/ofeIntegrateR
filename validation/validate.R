##############################################################################
## validation/validate.R
##
## Accuracy checks for ofeIntegrateR against independent references:
##   * ASReml-R, fitting the same model to the same data (the REML engine,
##     means, SEDs and Wald tests);
##   * gstat called directly with the true variogram (the kriging functions);
##   * a known truth from simulation (bias, interval coverage, type I error,
##     variogram recovery, defect detection);
##   * structural invariants (gridding, designs, sampling, letter displays).
##
## Unit tests check that functions run and return the right shape. This
## script checks that their answers are right, so it is slower and is not run
## by R CMD check. Run from the package root:
##
##   Rscript validation/validate.R [n_cores]
##
## Writes validation/results/checks.csv (one row per check, with pass/fail)
## and validation/results/montecarlo.csv.
##############################################################################

suppressWarnings(suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(gstat); library(sp)
}))
has_asreml <- requireNamespace("asreml", quietly = TRUE)
has_sommer <- requireNamespace("sommer", quietly = TRUE)
if (has_asreml) {
  suppressMessages(library(asreml))
  asreml::asreml.options(trace = FALSE)
}
args   <- commandArgs(trailingOnly = TRUE)
N_CORE <- if (length(args)) as.integer(args[1]) else 1L
dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)

checks <- list()
record <- function(area, check, value, reference, tol, pass = NULL, note = "") {
  if (is.null(pass)) pass <- isTRUE(abs(value - reference) <= tol)
  checks[[length(checks) + 1L]] <<- data.frame(
    area = area, check = check, value = signif(value, 5),
    reference = signif(reference, 5), tolerance = tol, pass = pass,
    note = note, stringsAsFactors = FALSE)
  cat(sprintf("[%s] %-14s %-58s %10.5g vs %10.5g\n",
              if (pass) "PASS" else "FAIL", area, check, value, reference))
}
relerr <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-8))
par_map <- function(X, FUN, ...) {
  if (N_CORE > 1) parallel::mclapply(X, FUN, ..., mc.cores = N_CORE)
  else lapply(X, FUN, ...)
}

## ===========================================================================
## 1. fit_ofe against asreml: same model, same data
## ===========================================================================
cat("\n== 1. REML engine: fit_ofe vs asreml ==\n")
if (has_asreml) {
  for (s in 1:3) {
    sim <- simulate_ofe_trial(n_row = 30, n_col = 18, seed = s)
    g <- sim$grid[order(sim$grid$row, sim$grid$col), ]
    g$rowf <- factor(g$row); g$colf <- factor(g$col)

    fo <- fit_ofe(dense_response ~ treat, data = g,
                  residual = ~ ar1(row):ar1(col))
    fa <- asreml::asreml(dense_response ~ treat,
                         residual = ~ ar1(rowf):ar1(colf), data = g)
    vo <- setNames(fo$varcomp$estimate, fo$varcomp$component)
    va <- summary(fa)$varcomp
    ca <- summary(fa, coef = TRUE)$coef.fixed

    record("fit_ofe", sprintf("AR1xAR1 seed %d: row correlation", s),
           vo[["R!row!cor"]], va["rowf:colf!rowf!cor", "component"], 1e-3)
    record("fit_ofe", sprintf("AR1xAR1 seed %d: col correlation", s),
           vo[["R!col!cor"]], va["rowf:colf!colf!cor", "component"], 1e-3)
    record("fit_ofe", sprintf("AR1xAR1 seed %d: residual variance (rel.)", s),
           relerr(vo[["sigma2"]], va["rowf:colf!R", "component"]), 0, 1e-3)
    record("fit_ofe", sprintf("AR1xAR1 seed %d: treatment effects (max rel.)", s),
           relerr(fo$coefficients[c("treatB", "treatC")],
                  ca[c("treat_B", "treat_C"), "solution"]), 0, 1e-3)
    record("fit_ofe", sprintf("AR1xAR1 seed %d: treatment SEs (max rel.)", s),
           relerr(sqrt(diag(fo$vcov))[c("treatB", "treatC")],
                  ca[c("treat_B", "treat_C"), "std error"]), 0, 1e-3)
  }

  # A random effect alongside the spatial residual: replicate blocks
  p <- simulate_paddock(n_row = 30, n_col = 24, treatments = c("A", "B", "C"),
                        n_rep = 4, seed = 7)
  p <- p[order(p$row, p$col), ]
  p$rowf <- factor(p$row); p$colf <- factor(p$col); p$rep <- factor(p$rep)
  fo <- fit_ofe(yield ~ treat, random = ~ rep, data = p,
                residual = ~ ar1(row):ar1(col))
  fa <- asreml::asreml(yield ~ treat, random = ~ rep,
                       residual = ~ ar1(rowf):ar1(colf), data = p)
  ca <- summary(fa, coef = TRUE)$coef.fixed
  record("fit_ofe", "random rep + AR1xAR1: treatment effects (max rel.)",
         relerr(fo$coefficients[c("treatB", "treatC")],
                ca[c("treat_B", "treat_C"), "solution"]), 0, 2e-3)
  record("fit_ofe", "random rep + AR1xAR1: treatment SEs (max rel.)",
         relerr(sqrt(diag(fo$vcov))[c("treatB", "treatC")],
                ca[c("treat_B", "treat_C"), "std error"]), 0, 2e-3)
  # asreml drops constants from its log-likelihood, so compare the change in
  # REML log-likelihood between two models with the same fixed effects
  fo0 <- fit_ofe(yield ~ treat, random = ~ rep, data = p,
                 residual = ~ id(row):id(col))
  fa0 <- asreml::asreml(yield ~ treat, random = ~ rep, data = p)
  record("fit_ofe", "random rep: REML logLik gain AR1xAR1 over id, vs asreml",
         fo$loglik - fo0$loglik, fa$loglik - fa0$loglik, 0.01)

  # Separate residual sections per pseudo-environment (dsum), via adaptive_residual()
  p2 <- simulate_paddock(n_row = 30, n_col = 24, n_zones = 2,
                         treatments = c("A", "B", "C"), n_rep = 4, seed = 11)
  p2 <- p2[order(p2$row, p2$col), ]
  rf <- adaptive_residual(p2, zone = "zone")
  fo <- fit_ofe(yield ~ zone + zone:treat, data = p2, residual = rf)
  p2$rowf <- factor(p2$row); p2$colf <- factor(p2$col)
  fa <- asreml::asreml(yield ~ zone + zone:treat,
                       residual = ~ dsum(~ ar1(rowf):ar1(colf) | zone), data = p2)
  va <- summary(fa)$varcomp
  vo <- fo$varcomp
  record("fit_ofe", "dsum by zone: section variances (max rel.)",
         relerr(sort(vo$variance[vo$type == "var"]),
                sort(va$component[grepl("!R$", rownames(va))])), 0, 5e-3)
  bo <- fo$coefficients[grepl("treat", names(fo$coefficients))]
  ba <- summary(fa, coef = TRUE)$coef.fixed
  ba <- ba[grepl("treat", rownames(ba)) & ba[, "solution"] != 0, "solution"]
  record("fit_ofe", "dsum by zone: zone:treat effects (max rel.)",
         relerr(sort(bo), sort(ba)), 0, 5e-3)

  # Independent residual: must equal ordinary least squares
  sim <- simulate_ofe_trial(n_row = 20, n_col = 12, seed = 3)
  fo <- fit_ofe(dense_response ~ treat, data = sim$grid,
                residual = ~ id(row):id(col))
  fl <- lm(dense_response ~ treat, data = sim$grid)
  record("fit_ofe", "id residual: coefficients equal lm (max abs.)",
         max(abs(fo$coefficients - coef(fl))), 0, 1e-8)
  record("fit_ofe", "id residual: SEs equal lm (max abs.)",
         max(abs(sqrt(diag(fo$vcov)) - sqrt(diag(vcov(fl))))), 0, 1e-8)
}

## ===========================================================================
## 2. Means, SEDs, Wald tests and letters
## ===========================================================================
cat("\n== 2. Means, SEDs, Wald tests, letters ==\n")
if (has_asreml) {
  sim <- simulate_ofe_trial(n_row = 30, n_col = 18, n_treat = 4,
                            treat_effects = c(0, 0.3, 0.6, 0.65), seed = 5)
  g <- sim$grid[order(sim$grid$row, sim$grid$col), ]
  g$rowf <- factor(g$row); g$colf <- factor(g$col)
  fo <- fit_ofe(dense_response ~ treat, data = g, residual = ~ ar1(row):ar1(col))
  fa <- asreml::asreml(dense_response ~ treat,
                       residual = ~ ar1(rowf):ar1(colf), data = g)
  mo <- ofe_means(fo, "treat")
  pa <- predict(fa, classify = "treat", sed = TRUE)
  record("ofe_means", "predicted means (max abs.)",
         max(abs(mo$estimate - pa$pvals$predicted.value)), 0, 1e-3)
  pw <- ofe_means(fo, "treat", pairwise = TRUE)
  sed_a <- as.matrix(pa$sed)
  lv <- levels(g$treat)
  sed_ref <- mapply(function(a, b) sed_a[match(a, lv), match(b, lv)],
                    pw$level1, pw$level2)
  record("ofe_means", "pairwise SEDs vs asreml predict (max rel.)",
         relerr(pw$std.error, sed_ref), 0, 2e-3)
  wo <- wald_tests(fo)
  wa <- asreml::wald(fa, denDF = "none")
  record("wald_tests", "treatment Wald statistic vs asreml (rel.)",
         relerr(wo$wald[wo$term == "treat"], wa["treat", "Wald statistic"]), 0, 5e-3)

  # Letters must agree with the pairwise tests they summarise: two means
  # share a letter exactly when their comparison is not significant
  lsd <- ofe_lsd(fo, "treat")
  cmp <- attr(lsd, "comparisons")
  letters_of <- setNames(lsd$group, lsd$treat)
  share <- function(a, b) length(intersect(strsplit(letters_of[[a]], "")[[1]],
                                           strsplit(letters_of[[b]], "")[[1]])) > 0
  agree <- mapply(function(a, b, pv) share(a, b) == (pv >= 0.05),
                  as.character(pw$level1), as.character(pw$level2), pw$p.value)
  record("ofe_lsd", "letters agree with pairwise p-values (share of pairs)",
         mean(agree), 1, 0)
}

## ===========================================================================
## 3. Operating characteristics against a known truth (Monte Carlo)
## ===========================================================================
cat("\n== 3. Monte Carlo: bias, coverage, type I error ==\n")
mc_one <- function(i, effect) {
  p <- simulate_paddock(n_row = 30, n_col = 24, treatments = c("A", "B", "C"),
                        n_rep = 4, treat_effect = c(0, effect, 2 * effect),
                        seed = 1000 + i)
  f <- tryCatch(fit_ofe(yield ~ treat, random = ~ rep,
                        data = transform(p, rep = factor(rep)),
                        residual = ~ ar1(row):ar1(col)),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  pw <- ofe_means(f, "treat", pairwise = TRUE)
  ca <- pw[pw$level1 == "A" & pw$level2 == "C", ]
  wt <- wald_tests(f)
  p_treat <- wt$p.value[wt$term == "treat"]
  data.frame(effect = effect, rep = i, est = ca$estimate, se = ca$std.error,
             truth = 2 * effect, p_contrast = ca$p.value, p_treat = p_treat)
}
N_MC <- 200L
mc <- do.call(rbind, c(par_map(seq_len(N_MC), mc_one, effect = 0),
                       par_map(seq_len(N_MC), mc_one, effect = 0.3)))
write.csv(mc, "validation/results/montecarlo.csv", row.names = FALSE)
m0 <- mc[mc$effect == 0, ]; m1 <- mc[mc$effect == 0.3, ]
record("fit_ofe MC", "null: type I error of treatment Wald test",
       mean(m0$p_treat < 0.05, na.rm = TRUE), 0.05, 0.03,
       note = sprintf("%d replicates; truth is an exponential field + noise, fitted as AR1xAR1", nrow(m0)))
record("fit_ofe MC", "C-A contrast: bias / true effect",
       mean(m1$est - m1$truth) / 0.6, 0, 0.05)
record("fit_ofe MC", "C-A contrast: 95% CI coverage",
       mean(abs(m1$est - m1$truth) <= qnorm(0.975) * m1$se), 0.95, 0.03)
record("fit_ofe MC", "C-A contrast: mean SE / empirical SD",
       sqrt(mean(m1$se^2)) / sd(m1$est), 1, 0.1)

## ===========================================================================
## 4. Geostatistics: variogram, kriging, cross-validation, sample interval
## ===========================================================================
cat("\n== 4. Variogram, kriging, sample interval ==\n")
# Truth: Exp model, nugget 0.2, partial sill 1, range parameter 10 (practical 30)
true_vgm <- gstat::vgm(psill = 1, model = "Exp", range = 10, nugget = 0.2)
field_grid <- expand.grid(col = 1:60, row = 1:60)
sim_field <- function(seed) {
  set.seed(seed)
  z <- gstat::gstat(formula = z ~ 1, locations = ~ col + row, dummy = TRUE,
                    beta = 0, model = gstat::vgm(psill = 1, model = "Exp",
                                                 range = 10), nmax = 40)
  f <- predict(z, newdata = field_grid, nsim = 1, debug.level = 0)
  data.frame(col = field_grid$col, row = field_grid$row, truth = f$sim1)
}
vg_one <- function(s) {
  f <- sim_field(s)
  set.seed(s + 1)
  pts <- f[sample(nrow(f), 250), ]
  pts$obs <- pts$truth + rnorm(nrow(pts), 0, sqrt(0.2))
  v <- tryCatch(ofe_variogram(pts, value = "obs", x = "col", y = "row",
                              model = "exponential"), error = function(e) NULL)
  if (is.null(v)) return(NULL)
  data.frame(seed = s, nugget = v$nugget, psill = v$psill,
             practical_range = v$practical_range)
}
vg <- do.call(rbind, par_map(1:30, vg_one))
record("ofe_variogram", "median practical range (true 30)",
       median(vg$practical_range), 30, 30 * 0.3)
record("ofe_variogram", "median nugget (true 0.2)", median(vg$nugget), 0.2, 0.1)
record("ofe_variogram", "median total sill (true 1.2)",
       median(vg$nugget + vg$psill), 1.2, 1.2 * 0.3,
       note = "one 60x60 field per replicate; sill estimates vary with the realised field")

# Kriging: package vs oracle kriging with the true variogram
kr_one <- function(s) {
  f <- sim_field(100 + s)
  set.seed(s)
  idx <- sample(nrow(f), 120)
  pts <- f[idx, ]; pts$obs <- pts$truth + rnorm(nrow(pts), 0, sqrt(0.2))
  target <- f[-idx, ]
  k <- tryCatch(krige_point_samples(pts, target, value = "obs",
                                    coords = c("col", "row")),
                error = function(e) NULL)
  if (is.null(k)) return(NULL)
  o <- gstat::krige(obs ~ 1, locations = ~ col + row, data = pts,
                    newdata = target, model = true_vgm, nmax = 30,
                    debug.level = 0)
  # The prediction target is the observable (truth + measurement error)
  y_obs <- target$truth + rnorm(nrow(target), 0, sqrt(0.2))
  z <- (y_obs - k$obs_kriged) / sqrt(k$obs_kriged_var)
  data.frame(seed = s,
             rmse_pkg = sqrt(mean((k$obs_kriged - target$truth)^2)),
             rmse_oracle = sqrt(mean((o$var1.pred - target$truth)^2)),
             cover95 = mean(abs(z) <= 1.96), z_sd = sd(z))
}
kr <- do.call(rbind, par_map(1:30, kr_one))
record("krige_point_samples", "RMSE / oracle-kriging RMSE (median)",
       median(kr$rmse_pkg / kr$rmse_oracle), 1, 0.1)
record("krige_point_samples", "95% prediction-interval coverage (median)",
       median(kr$cover95), 0.95, 0.04)
record("krige_point_samples", "SD of standardised prediction errors (median)",
       median(kr$z_sd), 1, 0.15)

# Cross-validation matches gstat::krige.cv with the same variogram
f <- sim_field(7); set.seed(7)
pts <- f[sample(nrow(f), 80), ]; pts$obs <- pts$truth + rnorm(80, 0, sqrt(0.2))
cv <- cv_krige_surface(pts, value = "obs", coords = c("col", "row"))
vm <- attr(cv, "variogram")
if (inherits(vm, "ofe_variogram")) vm <- vm$fits[[1]] %||% vm$model
if (inherits(vm, "variogramModel")) {
  ref <- gstat::krige.cv(obs ~ 1, locations = ~ col + row, data = pts,
                         model = vm, nmax = 30, verbose = FALSE)
  record("cv_krige_surface", "RMSE equals gstat::krige.cv",
         cv$rmse, sqrt(mean(ref$residual^2)), 1e-8)
} else {
  record("cv_krige_surface", "variogram attached", 0, 1, 0,
         note = "attr(, 'variogram') is not a gstat model; compare by hand")
}

# Sample interval: reproduce the design calculation directly with gstat
si <- kriging_sample_interval(nugget = 0.2, psill = 1, range = 30,
                              target_kse = 0.7, area_ha = 10)
d <- si$interval
corner <- data.frame(x = c(0, d, 0, d), y = c(0, 0, d, d), z = 0)
centre <- data.frame(x = d / 2, y = d / 2)
m_prac <- gstat::vgm(psill = 1, model = "Exp", range = 30 / 3, nugget = 0.2)
kv <- gstat::krige(z ~ 1, locations = ~ x + y, data = corner,
                   newdata = centre, model = m_prac, debug.level = 0)$var1.var
record("kriging_sample_interval", "relative kriging SE at chosen interval",
       si$rel_kse, sqrt(kv / 1.2), 0.01,
       note = "range passed as practical range (30); gstat Exp range parameter = 30/3")
record("kriging_sample_interval", "floor equals sqrt(nugget / sill)",
       si$kse_floor, sqrt(0.2 / 1.2), 1e-6)
record("kriging_sample_interval", "achieved rel. SE does not exceed target",
       si$rel_kse, 0.7, 0, pass = si$rel_kse <= 0.7 + 1e-9)

## ===========================================================================
## 5. Gridding and cleaning of yield-monitor data
## ===========================================================================
cat("\n== 5. Gridding and cleaning ==\n")
ym <- simulate_yield_monitor(n_point_samples = 20, seed = 3)
g <- grid_dense_layer(ym$cloud, response = "yield", treat = "treat",
                      cell_size = 9)
record("grid_dense_layer", "observations conserved (sum n_obs = nrow cloud)",
       sum(g$n_obs, na.rm = TRUE), nrow(ym$cloud), 0)
record("grid_dense_layer", "n_obs-weighted mean of cell means = cloud mean",
       sum(g$yield * g$n_obs, na.rm = TRUE) / sum(g$n_obs, na.rm = TRUE),
       mean(ym$cloud$yield), 1e-8)
record("grid_dense_layer", "complete lattice (rows x cols = cells)",
       nrow(g), length(unique(g$row)) * length(unique(g$col)), 0)

# Settings the function documents for this simulator: pass_trim = 12 (the
# simulator damages 12 m at each pass end), a speed column, local_mad = 3.
# min_yield cannot be used here because the simulated yield is on a
# standardised scale where 0 is an ordinary value, so zero-reading dropouts
# are only partly caught; on data in t/ha, min_yield = 0 removes them.
cl_rate <- t(sapply(1:10, function(s) {
  yd <- simulate_yield_monitor(defects = TRUE, seed = s)
  d <- yd$cloud$defect; d[is.na(d)] <- ""
  cl <- clean_yield_monitor(yd$cloud, yield = "yield", pass = "pass",
                            speed = "speed", pass_trim = 12, local_mad = 3,
                            drop = FALSE)
  flag <- !(cl$.keep %in% TRUE)
  c(ends = mean(flag[d %in% c("pass_start", "pass_end")]),
    outlier = mean(flag[d == "outlier"]), good = mean(flag[d == ""]))
}))
record("clean_yield_monitor", "recall: pass-start and pass-end defects (10 fields)",
       mean(cl_rate[, "ends"]), 1, 0, pass = mean(cl_rate[, "ends"]) >= 0.9)
record("clean_yield_monitor", "sound points removed (10 fields)",
       mean(cl_rate[, "good"]), 0, 0, pass = mean(cl_rate[, "good"]) <= 0.02)
record("clean_yield_monitor", "recall: outliers, without min_yield (reported)",
       mean(cl_rate[, "outlier"]), 1, 1, pass = TRUE,
       note = "documented: ~39% without min_yield, ~94% with min_yield = 0 on t/ha data")

## ===========================================================================
## 6. Integration fits: engines must agree
## ===========================================================================
cat("\n== 6. Integration fits ==\n")
sim <- simulate_ofe_trial(n_row = 30, n_col = 18, n_point_samples = 40, seed = 9)
kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
kr <- kr[order(kr$row, kr$col), ]
if (has_asreml) {
  a <- extract_fixed_effects(fit_integrated_kriged(kr, "dense_response", "treat",
                                                   "point_obs_kriged", engine = "asreml"))
  o <- extract_fixed_effects(fit_integrated_kriged(kr, "dense_response", "treat",
                                                   "point_obs_kriged", engine = "ofe"))
  get <- function(fx, pat) fx$estimate[grepl(pat, fx$term)][1]
  record("fit_integrated_kriged", "treatment C: ofe engine vs asreml engine (rel.)",
         relerr(get(o, "treat_?C$"), get(a, "treat_?C$")), 0, 2e-3)
  record("fit_integrated_kriged", "covariate slope: ofe vs asreml (rel.)",
         relerr(get(o, "point_obs_kriged"), get(a, "point_obs_kriged")), 0, 2e-3)
  ci <- compare_integration(kr, "dense_response", "treat", "point_obs_kriged",
                            engine = "ofe")
  record("compare_integration", "returns baseline and integrated estimates",
         ncol(ci), 1, 0, pass = nrow(ci) > 0 && !is.null(attr(ci, "models")))
}
if (has_asreml && has_sommer) {
  # sommer only accepts two complete layers (it imputes missing responses)
  both <- sim$grid[, c("row", "col", "point_true")]
  ja <- extract_fixed_effects(fit_integrated_joint(
    sim$grid, both, "dense_response", "point_true", residual = "id",
    min_overlap = 0))
  js <- extract_fixed_effects(fit_integrated_joint(
    sim$grid, both, "dense_response", "point_true", engine = "sommer",
    min_overlap = 0))
  record("fit_integrated_joint", "complete layers: sommer vs asreml, treatment C (rel.)",
         relerr(js$estimate[js$term == "dense_response:treatC"],
                ja$estimate[ja$term == "at(trait, dense_response):treat_C"]), 0, 1e-3)
  sp_err <- tryCatch({fit_integrated_joint(sim$grid, sim$point_samples,
                        "dense_response", "point_obs", engine = "sommer"); FALSE},
                     error = function(e) TRUE)
  record("fit_integrated_joint", "sommer refuses a sparse layer", as.numeric(sp_err), 1, 0)
}
if (has_asreml) {
  # Sparse response recovered against a known truth (the M4 follow-up case):
  # hand cuts at 4 cells per plot, dense NDVI that responds to treatment
  jt_one <- function(i) {
    p <- simulate_paddock(n_row = 24, n_col = 18, treatments = c("N0", "N60", "N120"),
                          n_rep = 3, treat_effect = c(0, 0.3, 0.5), seed = 500 + i)
    set.seed(500 + i)
    p$ndvi <- 0.55 + 0.06 * as.numeric(scale(p$yield_potential)) +
      c(N0 = 0, N60 = 0.03, N120 = 0.05)[as.character(p$treat)] +
      rnorm(nrow(p), 0, 0.02)
    cuts <- do.call(rbind, lapply(split(p, p$plot), function(d) d[sample(nrow(d), 4), ]))
    f <- tryCatch(suppressWarnings(fit_integrated_joint(
      p[, c("row", "col", "treat", "rep", "ndvi")],
      data.frame(row = cuts$row, col = cuts$col, cut_yield = cuts$yield),
      "ndvi", "cut_yield", block = "rep")), error = function(e) NULL)
    if (is.null(f) || !isTRUE(f$converge)) return(NULL)
    fx <- extract_fixed_effects(f)
    r <- fx[fx$term == "at(trait, cut_yield):treat_N120", ]
    u <- summary(lm(yield ~ factor(rep) + treat, data = cuts))$coefficients["treatN120", ]
    data.frame(est = r$estimate, se = r$se, uni_est = u[["Estimate"]])
  }
  jt <- do.call(rbind, par_map(1:60, jt_one))
  record("fit_integrated_joint", "sparse response: bias of N120-N0 / effect",
         mean(jt$est - 0.5) / 0.5, 0, 0.1, note = sprintf("%d converged fits", nrow(jt)))
  record("fit_integrated_joint", "sparse response: 95% CI coverage",
         mean(abs(jt$est - 0.5) <= 1.96 * jt$se), 0.95, 0.06)
  record("fit_integrated_joint", "sparse response: empirical SD joint / cuts alone",
         sd(jt$est) / sd(jt$uni_est), 1, 0, pass = sd(jt$est) < sd(jt$uni_est),
         note = "below 1 means the dense layer made the sparse response more precise")
}

## ===========================================================================
## 7. Designs, sampling, zoning
## ===========================================================================
cat("\n== 7. Designs, sampling, zoning ==\n")
des <- as.data.frame(make_trial_design(c("N0", "N60", "N120"), n_rep = 5,
                                       seed = 1))
per_rep <- tapply(des$treat, des$rep, function(t) length(unique(t)))
record("make_trial_design", "every replicate holds every treatment",
       min(per_rep), 3, 0)

sim <- simulate_ofe_trial(n_row = 30, n_col = 18, seed = 2)
for (d in c("random", "grid", "stratified", "nested")) {
  set.seed(1)
  pp <- tryCatch(place_point_samples(sim$grid, n = 24, design = d),
                 error = function(e) NULL)
  ok <- !is.null(pp) && all(paste(pp$row, pp$col) %in% paste(sim$grid$row, sim$grid$col))
  record("place_point_samples", sprintf("%s: points lie on the grid", d),
         as.numeric(ok), 1, 0)
  if (!is.null(pp)) record("place_point_samples", sprintf("%s: n returned", d),
                           nrow(pp), 24, 0, pass = nrow(pp) <= 24 && nrow(pp) >= 20,
                           note = "grid and nested designs may round the count")
}

# Zones: partition_paddock on covariates built with a known zone structure
set.seed(3)
zg <- expand.grid(col = 1:30, row = 1:30)
zg$zone_true <- factor(ifelse(zg$col <= 12, 1, ifelse(zg$row <= 15, 2, 3)))
mu <- c(0, 2, 4)[zg$zone_true]
zg$elev <- mu + rnorm(nrow(zg), 0, 0.5)
zg$clay <- c(1, -1, 1)[zg$zone_true] * 1.5 + rnorm(nrow(zg), 0, 0.5)
zp <- partition_paddock(zg, covariates = c("elev", "clay"), k = 3,
                        row = "row", col = "col")
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  s <- sum(choose(tab, 2)); sa <- sum(choose(rowSums(tab), 2))
  sb <- sum(choose(colSums(tab), 2)); e <- sa * sb / choose(n, 2)
  (s - e) / ((sa + sb) / 2 - e)
}
record("partition_paddock", "adjusted Rand index vs known zones (k = 3)",
       ari(zp$zone, zg$zone_true), 1, 0.1, pass = ari(zp$zone, zg$zone_true) >= 0.9)

## ===========================================================================
res <- do.call(rbind, checks)
write.csv(res, "validation/results/checks.csv", row.names = FALSE)
cat(sprintf("\n%d checks: %d pass, %d fail\n", nrow(res), sum(res$pass),
            sum(!res$pass)))
if (any(!res$pass)) print(res[!res$pass, c("area", "check", "value", "reference")])
