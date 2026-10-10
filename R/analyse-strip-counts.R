#' Compare treatments using hand counts, one mean per strip
#'
#' Analyses quadrat counts (plants, weeds, diseased plants) from a strip trial
#' at the level of the experimental unit: the counts in each strip are pooled
#' into a strip mean density, and treatments are compared on the log scale with
#' replicate blocks, giving treatment ratios with confidence intervals.
#'
#' Quadrats within a strip share its seeder pass and its position, so they are
#' not independent replicates of the treatment. Treating them as replicates
#' gave "significant" differences in up to 70% of simulated trials with no
#' treatment effect (AAGI-CU-RD-OFE Milestone 5); the strip-mean analysis kept
#' this at the nominal 5%. With two treatments and blocks, the analysis here is
#' the paired t-test on log strip means.
#'
#' When `baseline` is given (counts at the same pegged quadrats before the
#' treatment was applied, e.g. weeds before spraying), the response is the log
#' of the after/before ratio for each strip. This removes most of the
#' patchiness and roughly doubled the power for weeds and disease in the
#' Milestone 5 simulations.
#'
#' @param data Data frame with one row per quadrat.
#' @param count Name of the column of counts.
#' @param strip Name of the column identifying the strip (plot).
#' @param treat Name of the treatment column.
#' @param block Optional name of the replicate (block) column.
#' @param area Name of a column of quadrat areas (m^2), or a single number for
#'   a common area. Densities are reported per m^2.
#' @param baseline Optional name of a column of counts made before treatment at
#'   the same quadrats.
#' @param level Confidence level for the ratios.
#'
#' @return An object of class `ofe_strip_counts`: a list with `strips` (one row
#'   per strip: counts, area, density and the analysed response), `anova` (the
#'   F test for treatment), `ratios` (each treatment against the first level,
#'   as a ratio of densities with a confidence interval) and `model` (the
#'   fitted [stats::lm()]).
#'
#' @seealso [place_count_stations()], [count_precision()]
#'
#' @examples
#' set.seed(1)
#' # 4 replicates x 2 treatments, 8 stations x 2 quadrats of 0.25 m2
#' strips <- data.frame(strip = 1:8, rep = rep(1:4, each = 2),
#'                      treat = rep(c("Control", "High rate"), 4))
#' strips$dens <- 120 * ifelse(strips$treat == "High rate", 1.2, 1) *
#'   exp(rnorm(8, 0, 0.1))
#' q <- strips[rep(1:8, each = 16), ]
#' q$count <- rpois(nrow(q), q$dens * 0.25)
#' fit <- analyse_strip_counts(q, count = "count", strip = "strip",
#'                             treat = "treat", block = "rep", area = 0.25)
#' fit
#'
#' @export
analyse_strip_counts <- function(data, count, strip, treat, block = NULL,
                                 area = 1, baseline = NULL, level = 0.95) {
  need <- c(count, strip, treat, block, baseline,
            if (is.character(area)) area)
  miss <- setdiff(need, names(data))
  if (length(miss)) {
    stop("Columns not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  a <- if (is.character(area)) data[[area]] else rep(area, nrow(data))
  if (any(!is.finite(a)) || any(a <= 0)) stop("Quadrat areas must be positive.", call. = FALSE)
  y <- data[[count]]
  if (any(y < 0, na.rm = TRUE)) stop("Counts must be non-negative.", call. = FALSE)
  ok <- !is.na(y) & !is.na(data[[strip]]) & !is.na(data[[treat]])
  if (!is.null(baseline)) ok <- ok & !is.na(data[[baseline]])
  d <- data[ok, , drop = FALSE]; y <- y[ok]; a <- a[ok]

  sid <- as.character(d[[strip]])
  strips <- data.frame(strip = unique(sid), stringsAsFactors = FALSE)
  idx <- match(strips$strip, sid)
  strips$treat <- d[[treat]][idx]
  if (!is.null(block)) strips$block <- d[[block]][idx]
  # every strip must carry a single treatment (and block)
  if (any(tapply(as.character(d[[treat]]), sid, function(v) length(unique(v))) > 1)) {
    stop("A strip has more than one treatment; check the `strip` column.", call. = FALSE)
  }
  strips$quadrats <- as.vector(table(factor(sid, levels = strips$strip)))
  strips$count <- as.vector(tapply(y, factor(sid, levels = strips$strip), sum))
  strips$area <- as.vector(tapply(a, factor(sid, levels = strips$strip), sum))
  strips$density <- strips$count / strips$area
  if (is.null(baseline)) {
    # 0.5 of a plant per strip keeps empty strips finite
    strips$response <- log((strips$count + 0.5) / strips$area)
  } else {
    b <- d[[baseline]]
    strips$baseline <- as.vector(tapply(b, factor(sid, levels = strips$strip), sum))
    strips$response <- log((strips$count + 0.5) / (strips$baseline + 0.5))
  }

  trt <- factor(strips$treat)
  if (nlevels(trt) < 2L) stop("At least two treatments are needed.", call. = FALSE)
  if (any(table(trt) < 2L) && is.null(block)) {
    warning("Some treatments have only one strip: no replication to test against.",
            call. = FALSE)
  }
  strips$treat <- trt
  form <- if (is.null(block)) {
    stats::as.formula("response ~ treat")
  } else {
    strips$block <- factor(strips$block)
    stats::as.formula("response ~ block + treat")
  }
  fit <- stats::lm(form, data = strips)
  if (fit$df.residual < 1L) {
    stop("No residual degrees of freedom: more replicate strips are needed.", call. = FALSE)
  }
  av <- stats::anova(fit)
  cf <- summary(fit)$coefficients
  rows <- paste0("treat", levels(trt)[-1])
  tq <- stats::qt(1 - (1 - level) / 2, fit$df.residual)
  ratios <- data.frame(treatment = levels(trt)[-1], reference = levels(trt)[1],
                       log_ratio = cf[rows, 1], se = cf[rows, 2],
                       ratio = exp(cf[rows, 1]),
                       lower = exp(cf[rows, 1] - tq * cf[rows, 2]),
                       upper = exp(cf[rows, 1] + tq * cf[rows, 2]),
                       p_value = cf[rows, 4], df = fit$df.residual)
  rownames(ratios) <- NULL
  out <- list(strips = strips, anova = av, ratios = ratios, model = fit,
              baseline = !is.null(baseline), level = level)
  class(out) <- "ofe_strip_counts"
  out
}

#' @export
print.ofe_strip_counts <- function(x, ...) {
  cat("Strip-level analysis of counts:", nrow(x$strips), "strips,",
      sum(x$strips$quadrats), "quadrats\n")
  if (x$baseline) cat("Response: log(after / before) per strip\n")
  else cat("Response: log strip mean density\n")
  ft <- x$anova["treat", ]
  cat(sprintf("Treatment: F = %.2f on %d and %d df, p = %.3g\n\n",
              ft[["F value"]], ft[["Df"]], x$model$df.residual, ft[["Pr(>F)"]]))
  r <- x$ratios
  cat(sprintf("Ratio to %s (%.0f%% CI):\n", r$reference[1], 100 * x$level))
  for (i in seq_len(nrow(r))) {
    cat(sprintf("  %-15s %.3f  (%.3f, %.3f)  p = %.3g\n", r$treatment[i], r$ratio[i],
                r$lower[i], r$upper[i], r$p_value[i]))
  }
  invisible(x)
}
