# The empirical variogram of a point-sampled variable, and a model fitted to it.
#
# Everything downstream of sampling depends on this: how far apart to put the
# cores, whether kriging the samples will reconstruct anything, and whether a
# spatial model is worth fitting at all. Until now the package kriged without
# ever showing the user the variogram it was kriging from.

#' Semivariance model functions, as a fraction of the partial sill
#' @keywords internal
#' @noRd
.vgm_shape <- function(h, a, model) {
  switch(model,
         exponential = 1 - exp(-h / a),
         spherical   = ifelse(h >= a, 1,
                              1.5 * (h / a) - 0.5 * (h / a)^3),
         gaussian    = 1 - exp(-(h / a)^2))
}

#' Practical range: the lag at which correlation has decayed to 0.05
#' @keywords internal
#' @noRd
.vgm_practical <- function(a, model) {
  switch(model, exponential = 3 * a, spherical = a, gaussian = sqrt(3) * a)
}

#' Fit one variogram model by constrained grid search over the range
#'
#' Given the range, nugget and partial sill enter linearly, so only the range
#' needs searching. Both are held non-negative -- an unconstrained fit will buy
#' an arbitrarily short range with a negative nugget -- and lags are weighted by
#' their pair count, so the noisy long lags do not set the answer.
#' @keywords internal
#' @noRd
.vgm_fit_one <- function(h, gamma, n_pair, model) {
  a_grid <- exp(seq(log(0.25 * min(h)), log(3 * max(h)), length.out = 200))
  out <- lapply(a_grid, function(a) {
    b <- .vgm_shape(h, a, model)
    X <- cbind(1, b)
    cf <- tryCatch(as.numeric(solve(crossprod(X * n_pair, X),
                                    crossprod(X * n_pair, gamma))),
                   error = function(e) c(NA_real_, NA_real_))
    if (any(!is.finite(cf)) || cf[2] < 0) return(list(sse = Inf))
    if (cf[1] < 0) {
      cf[1] <- 0
      cf[2] <- sum(n_pair * b * gamma) / sum(n_pair * b * b)
    }
    list(sse = sum(n_pair * (gamma - (cf[1] + cf[2] * b))^2),
         nugget = cf[1], psill = cf[2], range = a)
  })
  sse <- vapply(out, function(z) z$sse, numeric(1))
  if (all(!is.finite(sse))) return(NULL)
  best <- out[[which.min(sse)]]
  best$model <- model
  best
}

#' Empirical and fitted variogram of a point-sampled variable
#'
#' Bins every pair of samples by the distance between them, averages the
#' squared difference within each bin, and fits a variogram model to the
#' result. This is the summary that decides the rest of a sampling programme:
#' the range says how far one core speaks for, and the nugget says how much of
#' the variation no sampling density will ever resolve.
#'
#' # Reading the output
#'
#' The **nugget ratio** -- nugget over total sill -- is the number to look at
#' first, and [print.ofe_variogram()] classifies it on the usual convention
#' (Cambardella et al. 1994): below 0.25 is strong spatial dependence, 0.25 to
#' 0.75 moderate, above 0.75 weak. A weak-dependence variable is effectively
#' noise at the scale sampled, and no amount of kriging will make a useful
#' surface of it; that is a finding about the sampling design, and the honest
#' response is to sample closer together or give up on that layer, not to krige
#' it anyway.
#'
#' # Removing a trend first
#'
#' A variogram assumes the mean is constant. If it is not -- the treatment
#' raises yield on half the strips, the paddock slopes -- the variogram keeps
#' climbing and never reaches a sill, and the fitted range is inflated by the
#' trend rather than describing the correlation. Pass `trend` to remove that
#' first: `trend = ~ treat` fits and strips the treatment effect, and the
#' variogram is computed on what remains.
#'
#' @param data Data frame of point samples.
#' @param value Character; the column to compute the variogram of.
#' @param x,y Character; projected coordinate columns, in metres. Defaults to
#'   `"x"`/`"y"`, or `"x_centre"`/`"y_centre"` when those are present.
#' @param trend Optional one-sided formula of terms to remove before computing
#'   the variogram, e.g. `~ treat` or `~ treat + rep`.
#' @param n_bin Integer; number of distance bins.
#' @param cutoff Numeric; largest lag to use. Defaults to a third of the
#'   greatest distance between samples, the usual convention -- beyond that the
#'   bins hold few pairs and say more about the paddock's shape than its soil.
#' @param model Variogram model(s) to fit: any of `"exponential"`,
#'   `"spherical"`, `"gaussian"`. With more than one, the best weighted fit is
#'   chosen and the others are kept for comparison.
#'
#' @return An object of class `ofe_variogram`: a list with `empirical` (a data
#'   frame of `h`, `gamma`, `n_pair`), the fitted `model`, `nugget`, `psill`,
#'   `range`, `practical_range`, `sill`, `nugget_ratio`, `range_identified`
#'   (whether the practical range falls inside the distances actually sampled),
#'   the `fits` for every model tried, and the settings used. [plot.ofe_variogram()] draws it;
#'   [kriging_sample_interval()] takes it directly to turn it into a sampling
#'   interval and a core count.
#'
#' @seealso [kriging_sample_interval()], [krige_point_samples()],
#'   [cv_krige_surface()].
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 30, n_col = 20, n_point_samples = 60,
#'                           point_range = 6, seed = 1)
#' pts <- sim$point_samples
#' pts$x <- pts$col
#' pts$y <- pts$row
#' v <- ofe_variogram(pts, value = "point_obs")
#' v
#' plot(v)
#'
#' @export
ofe_variogram <- function(data, value, x = NULL, y = NULL, trend = NULL,
                          n_bin = 15L, cutoff = NULL,
                          model = c("exponential", "spherical", "gaussian")) {
  model <- match.arg(model, several.ok = TRUE)
  data <- as.data.frame(data)
  if (missing(value) || length(value) != 1L || !is.character(value)) {
    stop("`value` must name a single column of `data`.", call. = FALSE)
  }
  if (is.null(x)) x <- if ("x" %in% names(data)) "x" else
    if ("x_centre" %in% names(data)) "x_centre" else NULL
  if (is.null(y)) y <- if ("y" %in% names(data)) "y" else
    if ("y_centre" %in% names(data)) "y_centre" else NULL
  if (is.null(x) || is.null(y)) {
    stop("Coordinates not found: pass `x` and `y`.", call. = FALSE)
  }
  needed <- c(value, x, y, if (!is.null(trend)) all.vars(trend))
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  keep <- stats::complete.cases(data[, needed, drop = FALSE])
  d <- data[keep, , drop = FALSE]
  n <- nrow(d)
  if (n < 10L) {
    stop("Only ", n, " complete samples: too few for a variogram. Around 30 ",
         "is a working minimum, and 50 or more is better.", call. = FALSE)
  }

  v <- d[[value]]
  if (!is.null(trend)) {
    f <- stats::as.formula(paste(".ofe_v ~", paste(all.vars(trend),
                                                   collapse = " + ")))
    d$.ofe_v <- v
    v <- stats::residuals(stats::lm(f, data = d))
  } else {
    v <- v - mean(v)
  }

  xy <- cbind(d[[x]], d[[y]])
  dist_all <- stats::dist(xy)
  if (is.null(cutoff)) cutoff <- max(dist_all) / 3
  g_all <- stats::dist(v)^2 / 2
  sel <- as.numeric(dist_all) > 0 & as.numeric(dist_all) <= cutoff
  if (sum(sel) < 30L) {
    stop("Only ", sum(sel), " sample pairs fall within the cutoff. Raise ",
         "`cutoff`, or sample more points.", call. = FALSE)
  }
  br <- seq(0, cutoff, length.out = as.integer(n_bin) + 1L)
  bin <- cut(as.numeric(dist_all)[sel], br, include.lowest = TRUE)
  gamma <- tapply(as.numeric(g_all)[sel], bin, mean)
  h_mid <- tapply(as.numeric(dist_all)[sel], bin, mean)
  n_pair <- tapply(as.numeric(g_all)[sel], bin, length)
  ok <- !is.na(gamma) & n_pair >= 2
  emp <- data.frame(h = as.numeric(h_mid[ok]), gamma = as.numeric(gamma[ok]),
                    n_pair = as.integer(n_pair[ok]))
  if (nrow(emp) < 3L) {
    stop("Fewer than three usable distance bins. Lower `n_bin`.", call. = FALSE)
  }

  fits <- lapply(model, function(m)
    .vgm_fit_one(emp$h, emp$gamma, emp$n_pair, m))
  names(fits) <- model
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (length(fits) == 0L) {
    stop("No variogram model could be fitted; the empirical variogram may be ",
         "flat. Check `plot()` of the empirical points.", call. = FALSE)
  }
  best <- fits[[which.min(vapply(fits, function(z) z$sse, numeric(1)))]]

  prac <- .vgm_practical(best$range, best$model)
  structure(list(
    empirical = emp, model = best$model, nugget = best$nugget,
    psill = best$psill, range = best$range,
    practical_range = prac,
    sill = best$nugget + best$psill,
    nugget_ratio = best$nugget / (best$nugget + best$psill),
    # A range longer than the largest lag observed is an extrapolation: the
    # variogram never reached its sill inside the data, so the fit is choosing
    # between curves that are indistinguishable over the distances sampled.
    # Reporting such a range without saying so is how a sampling plan ends up
    # built on a number the data never supported.
    range_identified = prac <= max(emp$h),
    sse = best$sse, fits = fits, n = n, value = value, cutoff = cutoff,
    max_lag = max(emp$h), trend = trend, coords = c(x, y),
    sd = stats::sd(d[[value]])
  ), class = "ofe_variogram")
}

#' Report a fitted variogram
#'
#' Prints the parameter table, classifies the spatial dependence from the
#' nugget ratio on the usual convention (Cambardella et al. 1994) and says what
#' that means for kriging the layer, and warns when the fitted range falls
#' outside the distances the samples actually cover.
#'
#' @param x An `ofe_variogram` object from [ofe_variogram()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly. Called for the printed report.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 30, n_col = 20, n_point_samples = 60,
#'                           point_range = 6, seed = 1)
#' pts <- sim$point_samples
#' pts$x <- pts$col
#' pts$y <- pts$row
#' print(ofe_variogram(pts, value = "point_obs"))
#'
#' @export
print.ofe_variogram <- function(x, ...) {
  cat("Variogram of `", x$value, "`  (", x$n, " samples, ",
      nrow(x$empirical), " lag bins)\n", sep = "")
  if (!is.null(x$trend)) {
    cat("Trend removed first: ", deparse(x$trend), "\n", sep = "")
  }
  cat("\nFitted model: ", x$model, "\n", sep = "")
  tab <- data.frame(
    value = c(x$nugget, x$psill, x$sill, x$range, x$practical_range,
              x$nugget_ratio),
    row.names = c("nugget (c0)", "partial sill (c1)", "total sill",
                  "range parameter", "practical range", "nugget ratio"))
  print(format(tab, digits = 4), quote = FALSE)

  if (!isTRUE(x$range_identified)) {
    cat("\n")
    cat(strwrap(paste0(
      "The fitted practical range (", format(x$practical_range, digits = 3),
      ") is longer than the largest lag the samples cover (",
      format(x$max_lag, digits = 3), "). The variogram has not reached its ",
      "sill inside the data, so the range is an extrapolation: curves that ",
      "differ well beyond the sampled distances fit these points equally. ",
      "Treat it as a lower bound, and sample over a longer transect before ",
      "using it to set a sampling interval."), width = 76, prefix = "  "),
      sep = "\n")
  }

  cls <- if (x$nugget_ratio < 0.25) "strong" else
    if (x$nugget_ratio <= 0.75) "moderate" else "weak"
  cat("\nSpatial dependence: ", cls, " (nugget ratio ",
      format(x$nugget_ratio, digits = 2), ")\n", sep = "")
  msg <- switch(
    cls,
    strong = paste0("Most of the variation is spatially structured; kriging ",
                    "this layer should reconstruct it well."),
    moderate = paste0("Part of the variation is spatially structured. Kriging ",
                      "will smooth, and the surface will be less variable ",
                      "than the truth."),
    weak = paste0("Most of the variation is unresolved at the distances ",
                  "sampled. Kriging this layer will return something close to ",
                  "its mean; sample closer together, or do not use it."))
  cat(strwrap(msg, width = 76, prefix = "  "), sep = "\n")
  if (length(x$fits) > 1L) {
    cat("\nModels compared (weighted SSE, lower is better):\n")
    ss <- vapply(x$fits, function(z) z$sse, numeric(1))
    print(format(data.frame(sse = ss[order(ss)]), digits = 4), quote = FALSE)
  }
  invisible(x)
}

#' Plot an empirical variogram and the model fitted to it
#'
#' Points are the empirical semivariances, sized by how many sample pairs each
#' bin rests on -- a bin built from six pairs should not be read as hard as one
#' built from six hundred, and drawing them the same size is the commonest way
#' a variogram plot misleads. The fitted model, the nugget and the practical
#' range are drawn over them.
#'
#' @param x An `ofe_variogram` from [ofe_variogram()].
#' @param show_all Logical; draw every fitted model rather than only the best.
#' @param main,xlab,ylab Plot labels.
#' @param ... Passed to [graphics::plot()].
#'
#' @return `x`, invisibly. Called for the plot.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 30, n_col = 20, n_point_samples = 60,
#'                           point_range = 6, seed = 1)
#' pts <- sim$point_samples
#' pts$x <- pts$col
#' pts$y <- pts$row
#' plot(ofe_variogram(pts, value = "point_obs"))
#'
#' @export
plot.ofe_variogram <- function(x, show_all = FALSE, main = NULL,
                               xlab = "distance", ylab = "semivariance", ...) {
  op <- .ofe_par(mar = c(4.6, 4.6, 3.2, 1.6))
  on.exit(graphics::par(op), add = TRUE)
  emp <- x$empirical
  if (is.null(main)) main <- paste0("Variogram of ", x$value)

  hh <- seq(0, max(emp$h) * 1.02, length.out = 200)
  curves <- if (show_all) x$fits else x$fits[x$model]
  ymax <- max(emp$gamma, vapply(curves, function(f)
    max(f$nugget + f$psill * .vgm_shape(hh, f$range, f$model)), numeric(1)))
  # Include the sill only when it is close enough to belong on the same axis;
  # stretching the axis to reach a sill the data never approached would shrink
  # the points that are the evidence.
  if (x$sill <= 1.3 * ymax) ymax <- max(ymax, x$sill)

  graphics::plot(emp$h, emp$gamma, type = "n", xlim = c(0, max(hh)),
                 ylim = c(0, ymax * 1.08), main = main, xlab = xlab,
                 ylab = ylab, xaxs = "i", yaxs = "i", ...)
  .ofe_grid(h = TRUE)

  # Reference lines only where they fall inside the data. Drawing a guide off
  # the edge of the panel tells the reader nothing; saying it is off the edge
  # tells them the range was not identified.
  if (x$sill <= ymax * 1.08) {
    graphics::abline(h = x$sill, col = .ofe_chrome[["axis"]], lty = 3)
  }
  if (isTRUE(x$range_identified)) {
    graphics::abline(v = x$practical_range, col = .ofe_chrome[["axis"]],
                     lty = 3)
  }

  cols <- ofe_palette("treatment", max(1L, length(curves)))
  for (i in seq_along(curves)) {
    f <- curves[[i]]
    graphics::lines(hh, f$nugget + f$psill * .vgm_shape(hh, f$range, f$model),
                    col = cols[i], lwd = 2)
  }

  # Area, not radius, proportional to pair count: radius would exaggerate.
  cexs <- 0.7 + 1.5 * sqrt(emp$n_pair / max(emp$n_pair))
  graphics::points(emp$h, emp$gamma, pch = 21, bg = .ofe_chrome[["surface"]],
                   col = .ofe_chrome[["ink2"]], cex = cexs, lwd = 1.4)

  # The model name goes in the subtitle, not beside the curve: at the right
  # edge the curve has flattened onto the sill and a label there lands on the
  # line and on the last points. In a small panel the full subtitle runs off
  # the edge, so shorten it rather than let it clip -- a truncated number is
  # worse than a missing one.
  flag <- if (isTRUE(x$range_identified)) "" else "  (extrapolated)"
  full <- sprintf("%s  |  nugget %.3g  |  sill %.3g  |  practical range %.3g%s",
                  if (show_all) "best fit, with alternatives" else x$model,
                  x$nugget, x$sill, x$practical_range, flag)
  short <- sprintf("%s  |  range %.3g%s", x$model, x$practical_range, flag)
  panel_w <- diff(graphics::par("usr")[1:2])
  sub <- if (graphics::strwidth(full, cex = 0.8) <= panel_w) full else
    if (graphics::strwidth(short, cex = 0.8) <= panel_w) short else x$model
  graphics::mtext(sub, side = 3, line = 0.25, adj = 0, cex = 0.8,
                  col = .ofe_chrome[["ink2"]])
  if (show_all && length(curves) > 1L) {
    graphics::legend("bottomright", legend = names(curves), col = cols, lwd = 2,
                     bty = "n", cex = 0.85, text.col = .ofe_chrome[["ink2"]])
  }
  invisible(x)
}
