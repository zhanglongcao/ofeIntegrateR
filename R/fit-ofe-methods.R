#' @export
print.ofe_fit <- function(x, ...) {
  cat("Spatial mixed model fitted by REML (ofeIntegrateR)\n\n")
  cat("Fixed:    ", deparse(x$fixed), "\n")
  if (!is.null(x$random)) cat("Random:   ", deparse(x$random), "\n")
  cat("Residual: ", deparse(x$residual), "\n")
  cat("\nObservations:", x$n, " Fixed parameters:", x$p,
      " REML logLik:", format(x$loglik, digits = 6), "\n")
  if (!x$converged) cat("*** the optimiser did not report convergence ***\n")
  invisible(x)
}

#' Summarise a fitted `ofe_fit`
#'
#' @param object An `ofe_fit` from [fit_ofe()].
#' @param ... Ignored.
#'
#' @return An object of class `summary.ofe_fit`, printed as a variance
#'   component table and a table of fixed-effect estimates. Standard errors on
#'   the variance parameters are asymptotic, from the observed information of
#'   the REML log-likelihood, and treat `sigma2` as independent of the
#'   correlation parameters; read them as a guide to identifiability rather
#'   than as exact inference.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 16, n_col = 8, seed = 1)
#' summary(fit_ofe(dense_response ~ treat, data = sim$grid))
#'
#' @export
summary.ofe_fit <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  tval <- object$coefficients / se
  coefs <- data.frame(
    term = names(object$coefficients),
    estimate = unname(object$coefficients),
    std.error = unname(se),
    statistic = unname(tval),
    p.value = 2 * stats::pt(abs(tval), df = object$n - object$p,
                            lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
  structure(list(fit = object, coefficients = coefs,
                 varcomp = object$varcomp),
            class = "summary.ofe_fit")
}

#' @export
print.summary.ofe_fit <- function(x, ...) {
  print(x$fit)
  cat("\nVariance parameters:\n")
  vc <- x$varcomp
  out <- data.frame(estimate = vc$estimate, std.error = vc$std.error,
                    variance = vc$variance, row.names = vc$component)
  print(format(out, digits = 4), quote = FALSE)
  cat("\nFixed effects:\n")
  cf <- x$coefficients
  out2 <- data.frame(estimate = cf$estimate, std.error = cf$std.error,
                     t = cf$statistic,
                     p = format.pval(cf$p.value, digits = 3, eps = 1e-4),
                     row.names = cf$term)
  print(format(out2, digits = 4), quote = FALSE)
  if (length(x$fit$aliased) > 0L) {
    cat("\nAliased (dropped from the model):",
        paste(x$fit$aliased, collapse = ", "), "\n")
  }
  if (x$fit$n_dropped > 0L) {
    cat("\n", x$fit$n_dropped,
        " row(s) dropped for missing values.\n", sep = "")
  }
  invisible(x)
}

#' @export
coef.ofe_fit <- function(object, ...) object$coefficients

#' @export
vcov.ofe_fit <- function(object, ...) object$vcov

#' @export
fitted.ofe_fit <- function(object, ...) object$fitted

#' @export
residuals.ofe_fit <- function(object, ...) object$residuals

#' @export
nobs.ofe_fit <- function(object, ...) object$n

#' @export
logLik.ofe_fit <- function(object, ...) {
  structure(object$loglik, df = nrow(object$varcomp), nobs = object$n,
            class = "logLik")
}

#' Wald tests of the fixed-effect terms
#'
#' The analogue of `wald.asreml()`: each term is tested by a Wald statistic on
#' all of its coefficients jointly, conditional on every other term in the
#' model. This is the test that answers "does treatment matter", as opposed to
#' the per-coefficient t-statistics in [summary.ofe_fit()].
#'
#' @param object An `ofe_fit` from [fit_ofe()].
#' @param ... Ignored.
#'
#' @return A data frame with one row per fixed term and columns `term`, `df`,
#'   `wald` (chi-square), `F` (`wald/df`), and `p.value`. The p-value uses an
#'   F reference distribution with `n - p` denominator degrees of freedom,
#'   which is approximate: it does not apply a Kenward-Roger correction, so it
#'   is mildly anti-conservative when the variance parameters are poorly
#'   determined.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 16, n_col = 8, seed = 1)
#' wald_tests(fit_ofe(dense_response ~ treat, data = sim$grid))
#'
#' @export
wald_tests <- function(object, ...) {
  UseMethod("wald_tests")
}

#' @rdname wald_tests
#' @export
wald_tests.ofe_fit <- function(object, ...) {
  labs <- attr(object$terms, "term.labels")
  asg <- object$assign
  beta <- object$coefficients
  V <- object$vcov
  den_df <- object$n - object$p

  rows <- lapply(seq_along(labs), function(j) {
    idx <- which(asg == j)
    if (length(idx) == 0L) return(NULL)
    b <- beta[idx]
    Vi <- V[idx, idx, drop = FALSE]
    inv <- tryCatch(solve(Vi), error = function(e) NULL)
    if (is.null(inv)) return(NULL)
    w <- drop(t(b) %*% inv %*% b)
    df <- length(idx)
    data.frame(term = labs[j], df = df, wald = w, F = w / df,
               p.value = stats::pf(w / df, df, den_df, lower.tail = FALSE),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  rownames(out) <- NULL
  out
}

# Build the contrast matrix L whose rows are the predicted cells: one row per
# level of `term`, or per level of `term` within each level of `by`. Shared by
# ofe_means() and ofe_lsd() so the means, the comparisons and the letters can
# never disagree with one another.
#' @keywords internal
#' @noRd
.ofe_cells <- function(object, term, by = NULL) {
  stopifnot(inherits(object, "ofe_fit"))
  dat <- object$data
  for (v in c(term, by)) {
    if (!v %in% names(dat)) {
      stop("`", v, "` is not a column of the fitted data.", call. = FALSE)
    }
  }
  lv <- object$xlevels[[term]]
  if (is.null(lv)) lv <- levels(factor(dat[[term]]))
  if (length(lv) < 2L) {
    stop("`", term, "` has fewer than two levels.", call. = FALSE)
  }

  if (is.null(by)) {
    groups <- list(seq_len(nrow(dat)))
    by_lv <- NA_character_
  } else {
    blv <- object$xlevels[[by]]
    if (is.null(blv)) blv <- levels(factor(dat[[by]]))
    bf <- factor(as.character(dat[[by]]), levels = blv)
    groups <- split(seq_len(nrow(dat)), bf)
    groups <- groups[lengths(groups) > 0L]
    by_lv <- names(groups)
    if (length(groups) == 0L) {
      stop("`", by, "` has no levels with data.", call. = FALSE)
    }
  }

  mt <- stats::delete.response(object$terms)
  keep_cols <- names(object$coefficients)
  rows <- list()
  lab_by <- character(0)
  lab_lv <- character(0)
  for (g in seq_along(groups)) {
    base <- dat[groups[[g]], , drop = FALSE]
    for (l in lv) {
      nd <- base
      nd[[term]] <- factor(rep(l, nrow(nd)), levels = lv)
      mm <- stats::model.matrix(mt, stats::model.frame(mt, nd,
                                                       xlev = object$xlevels))
      rows[[length(rows) + 1L]] <- colMeans(mm[, keep_cols, drop = FALSE])
      lab_by <- c(lab_by, if (is.null(by)) NA_character_ else by_lv[g])
      lab_lv <- c(lab_lv, l)
    }
  }
  L <- do.call(rbind, rows)
  colnames(L) <- keep_cols

  list(L = L, level = lab_lv, by = lab_by, levels = lv,
       by_levels = if (is.null(by)) NULL else by_lv,
       estimate = drop(L %*% object$coefficients),
       V = L %*% object$vcov %*% t(L),
       den_df = object$n - object$p)
}

#' Multiplicity adjustment for a set of pairwise comparisons
#' @keywords internal
#' @noRd
.ofe_adjust_p <- function(p, tstat, method, n_means, df) {
  switch(
    method,
    none  = p,
    tukey = stats::ptukey(abs(tstat) * sqrt(2), nmeans = n_means, df = df,
                          lower.tail = FALSE),
    sidak = pmin(1, 1 - (1 - p)^length(p)),
    stats::p.adjust(p, method = method)
  )
}

#' @keywords internal
#' @noRd
.ofe_adjust_methods <- c("none", "tukey", "sidak", "bonferroni", "holm",
                         "hochberg", "hommel", "BH", "BY", "fdr")

#' Pairwise comparisons within one predicted-cell group
#' @keywords internal
#' @noRd
.ofe_pairs <- function(est, V, lv, den_df, adjust) {
  cmb <- utils::combn(seq_along(lv), 2)
  d <- est[cmb[2, ]] - est[cmb[1, ]]
  sed <- sqrt(V[cbind(cmb[1, ], cmb[1, ])] + V[cbind(cmb[2, ], cmb[2, ])] -
                2 * V[cbind(cmb[1, ], cmb[2, ])])
  tv <- d / sed
  praw <- 2 * stats::pt(abs(tv), den_df, lower.tail = FALSE)
  data.frame(level1 = lv[cmb[1, ]], level2 = lv[cmb[2, ]],
             contrast = paste(lv[cmb[2, ]], "-", lv[cmb[1, ]]),
             estimate = unname(d), std.error = unname(sed),
             statistic = unname(tv),
             p.value = unname(.ofe_adjust_p(praw, tv, adjust, length(lv),
                                            den_df)),
             stringsAsFactors = FALSE)
}

#' Predicted means for a fixed term
#'
#' The analogue of `predict.asreml()`: the model-based mean of each level of a
#' factor, averaged over the other fixed effects as they actually occur in the
#' data (covariates are held at their observed mean). With `pairwise = TRUE`
#' it returns the differences between levels instead, each with the standard
#' error of the contrast -- which is the number a trial report quotes.
#'
#' @param object An `ofe_fit` from [fit_ofe()].
#' @param term Character; the name of a factor in the fixed model.
#' @param by Character or `NULL`; a second factor to compute the means within,
#'   one set per level of it. Use this for a pseudo-environment analysis --
#'   `ofe_means(fit, "treat", by = "zone")` gives the treatment means of each
#'   zone. Means are averaged over the rows belonging to that level only, and
#'   comparisons are made within a level, never across.
#' @param pairwise Logical; return all pairwise differences rather than the
#'   means themselves.
#' @param level Numeric; confidence level for the interval (default 0.95).
#' @param adjust Multiplicity adjustment applied to the pairwise p-values:
#'   `"none"` (the default, i.e. unprotected LSD comparisons), `"tukey"`,
#'   `"sidak"`, or any method taken by [stats::p.adjust()]. Ignored when
#'   `pairwise = FALSE`. With `by`, the adjustment counts the comparisons
#'   within a level, not across all of them.
#'
#' @return A data frame of predicted means (the `by` level where given, the
#'   `term` level, `estimate`, `std.error`, `lower`, `upper`) or, when
#'   `pairwise = TRUE`, of differences (`level1`, `level2`, `contrast`,
#'   `estimate`, `std.error`, `statistic`, `p.value`).
#'
#' @seealso [ofe_lsd()], which adds the a/b/c letters to the means.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 16, n_col = 8, seed = 1)
#' fit <- fit_ofe(dense_response ~ treat, data = sim$grid)
#' ofe_means(fit, "treat")
#' ofe_means(fit, "treat", pairwise = TRUE)
#' ofe_means(fit, "treat", pairwise = TRUE, adjust = "tukey")
#'
#' @export
ofe_means <- function(object, term, by = NULL, pairwise = FALSE, level = 0.95,
                      adjust = "none") {
  adjust <- match.arg(adjust, .ofe_adjust_methods)
  cells <- .ofe_cells(object, term, by)
  lv <- cells$levels
  den_df <- cells$den_df

  grp <- if (is.null(by)) rep(1L, length(cells$level)) else
    match(cells$by, cells$by_levels)

  if (!pairwise) {
    se <- sqrt(diag(cells$V))
    q <- stats::qt(1 - (1 - level) / 2, df = den_df)
    out <- data.frame(level = cells$level, estimate = unname(cells$estimate),
                      std.error = unname(se),
                      lower = unname(cells$estimate - q * se),
                      upper = unname(cells$estimate + q * se),
                      stringsAsFactors = FALSE)
    names(out)[1] <- term
    if (!is.null(by)) {
      out <- cbind(stats::setNames(data.frame(cells$by,
                                              stringsAsFactors = FALSE), by),
                   out)
    }
    return(out)
  }

  parts <- lapply(sort(unique(grp)), function(g) {
    i <- which(grp == g)
    tab <- .ofe_pairs(cells$estimate[i], cells$V[i, i, drop = FALSE],
                      cells$level[i], den_df, adjust)
    if (!is.null(by)) {
      tab <- cbind(stats::setNames(data.frame(rep(cells$by[i][1], nrow(tab)),
                                              stringsAsFactors = FALSE), by),
                   tab)
    }
    tab
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}
