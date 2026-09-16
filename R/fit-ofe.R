#' Control parameters for [fit_ofe()]
#'
#' @param maxit Integer; maximum iterations for each optimiser stage.
#' @param tol Numeric; convergence tolerance passed to [stats::optim()].
#' @param trace Logical; print the REML log-likelihood as it is optimised.
#' @param max_n Integer; refuse to fit more than this many observations. The
#'   engine is dense: it forms and factorises an `n` by `n` matrix at every
#'   iteration, so cost grows as `n^3`. Aggregate to a coarser lattice with
#'   [grid_dense_layer()] rather than raising this without thinking.
#' @param start Optional named numeric vector of starting values on the
#'   natural scale, named as in the `varcomp` table (correlations, ranges,
#'   relative variances). Useful when a fit has trouble converging.
#' @param se Logical; compute asymptotic standard errors for the variance
#'   parameters from the observed information. Costs one numerical Hessian.
#'
#' @return A list of control settings.
#' @export
ofe_control <- function(maxit = 500L, tol = 1e-8, trace = FALSE,
                        max_n = 2500L, start = NULL, se = TRUE) {
  list(maxit = as.integer(maxit), tol = tol, trace = isTRUE(trace),
       max_n = as.integer(max_n), start = start, se = isTRUE(se))
}

# Natural <-> unconstrained scale. Correlations live in (-1, 1), ranges and
# variance ratios on (0, Inf), and the optimiser sees only R^k.
#' @keywords internal
#' @noRd
.ofe_to_nat <- function(t, ptab, scales) {
  out <- numeric(length(t))
  for (i in seq_along(t)) {
    out[i] <- switch(ptab$type[i],
                     cor   = 0.9995 * tanh(t[i]),
                     range = scales[i] * exp(t[i]),
                     var   = exp(t[i]),
                     ratio = exp(t[i]))
  }
  out
}

#' @keywords internal
#' @noRd
.ofe_to_raw <- function(nat, ptab, scales) {
  out <- numeric(length(nat))
  for (i in seq_along(nat)) {
    out[i] <- switch(ptab$type[i],
                     cor   = atanh(pmax(-0.99, pmin(0.99, nat[i])) / 0.9995),
                     range = log(nat[i] / scales[i]),
                     var   = log(nat[i]),
                     ratio = log(nat[i]))
  }
  out
}

# Derivative of the natural parameter with respect to the unconstrained one,
# for the delta-method standard errors reported by summary().
#' @keywords internal
#' @noRd
.ofe_jacobian <- function(t, ptab, scales) {
  vapply(seq_along(t), function(i) {
    switch(ptab$type[i],
           cor   = 0.9995 / cosh(t[i])^2,
           range = scales[i] * exp(t[i]),
           var   = exp(t[i]),
           ratio = exp(t[i]))
  }, numeric(1))
}

#' Residual (restricted) log-likelihood at one parameter vector
#'
#' Returns `-Inf`-like penalties rather than errors when the covariance is not
#' positive definite, so the optimiser can walk back out of a bad region.
#' @keywords internal
#' @noRd
.ofe_eval <- function(theta, cache, full = FALSE) {
  n <- cache$n
  H <- .residual_matrix(cache$spec, theta, n)
  if (cache$n_random > 0L) {
    for (k in seq_len(cache$n_random)) {
      H <- H + theta[cache$ran_par[k]] * cache$ZZt[[k]]
    }
  }
  ch <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(ch)) return(if (full) NULL else 1e10)
  logdet_H <- 2 * sum(log(diag(ch)))

  Hiy <- backsolve(ch, forwardsolve(ch, cache$y, upper.tri = TRUE,
                                    transpose = TRUE))
  HiX <- backsolve(ch, forwardsolve(ch, cache$X, upper.tri = TRUE,
                                    transpose = TRUE))
  XtHiX <- crossprod(cache$X, HiX)
  chx <- tryCatch(chol(XtHiX), error = function(e) NULL)
  if (is.null(chx)) return(if (full) NULL else 1e10)
  logdet_XtHiX <- 2 * sum(log(diag(chx)))

  XtHiy <- crossprod(cache$X, Hiy)
  beta <- backsolve(chx, forwardsolve(chx, XtHiy, upper.tri = TRUE,
                                      transpose = TRUE))
  rss <- sum(cache$y * Hiy) - sum(XtHiy * beta)
  nu <- cache$n - cache$p
  if (!is.finite(rss) || rss <= 0) return(if (full) NULL else 1e10)
  sigma2 <- rss / nu

  # REML log-likelihood with sigma^2 profiled out.
  ll <- -0.5 * (logdet_H + logdet_XtHiX + nu * log(sigma2) + nu +
                  nu * log(2 * pi))
  if (!full) return(-ll)

  vcov_beta <- sigma2 * chol2inv(chx)
  list(loglik = ll, beta = drop(beta), vcov = vcov_beta, sigma2 = sigma2,
       H = H, chol_H = ch)
}

#' Fit a spatial mixed model without asreml
#'
#' A residual maximum likelihood (REML) fit of a linear mixed model with an
#' ASReml-style separable residual structure, written in base R. It exists so
#' that the analysis at the end of the OFE pipeline -- a treatment model with
#' an `ar1(row):ar1(col)` residual -- can be run by a collaborator who has no
#' ASReml-R licence, using the same formula spelling.
#'
#' # Residual structures
#'
#' `residual` accepts the asreml spellings, combined with `:` for a separable
#' (Kronecker) structure:
#'
#' \describe{
#'   \item{`id(f)`}{Independent; a bare variable name means the same thing.}
#'   \item{`ar1(f)`}{First-order autoregressive across the ordered, equally
#'     spaced levels of `f`. A factor uses its level order; a numeric column is
#'     read as a position on a lattice whose spacing is the commonest gap
#'     between its values, so a row missing from the data still counts as a
#'     lag rather than being closed up.}
#'   \item{`exp(x)`}{Exponential correlation in a numeric coordinate `x`,
#'     \eqn{\rho = \exp(-d/\phi)}. Use it when positions are irregular.}
#'   \item{`diag(f)`}{Heterogeneous variance across the levels of `f`, with the
#'     first level as the reference.}
#'   \item{`dsum(~ struct | s, levels = )`}{Independent sections: each named
#'     level of `s` gets its own copy of `struct` with its own parameters and
#'     its own variance. This is how a pseudo-environment analysis is written;
#'     [adaptive_residual()] builds the formula for you.}
#' }
#'
#' Unlike asreml, cells with a missing response are simply dropped: the
#' correlation is evaluated from the row and column positions of the
#' observations that remain, so an incomplete lattice needs no padding.
#'
#' # What it is not
#'
#' This is a compact, dense implementation intended for OFE-sized problems
#' (a few thousand lattice cells). It does not implement the sparse average
#' information algorithm, `us()`/`fa()` structures, or multi-trait models. For
#' those, use asreml where licensed; [fit_integrated_kriged()] will dispatch to
#' it with the same arguments.
#'
#' @param fixed Two-sided formula for the fixed effects, e.g.
#'   `yield ~ treat` or `yield ~ zone + zone:treat`.
#' @param data Data frame holding every variable used by `fixed`, `random` and
#'   `residual`. Rows with missing values in any of them are dropped.
#' @param random Optional one-sided formula of random effects, e.g. `~ rep` or
#'   `~ rep + rep:strip`. Each term contributes one variance component; all
#'   variables in a term are treated as factors.
#' @param residual One-sided formula giving the residual structure (see
#'   above). Defaults to `~ ar1(row):ar1(col)`.
#' @param control A list from [ofe_control()].
#'
#' @return An object of class `ofe_fit`: a list with elements `coefficients`,
#'   `vcov`, `sigma2`, `varcomp`, `loglik`, `fitted`, `residuals`, `n`, `p`,
#'   `converged`, and the model frame pieces needed by the methods. See
#'   [summary.ofe_fit()], [wald_tests()] and [ofe_means()].
#'
#' @seealso [adaptive_residual()] to build a `dsum()` formula from the realised
#'   geometry, [partition_pseudo_env()] to derive the sections in the first
#'   place, and [fit_integrated_kriged()] for the integration wrapper.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 16, n_col = 8, seed = 1)
#' fit <- fit_ofe(dense_response ~ treat, data = sim$grid,
#'                residual = ~ ar1(row):ar1(col))
#' summary(fit)
#' wald_tests(fit)
#' ofe_means(fit, "treat")
#'
#' @export
fit_ofe <- function(fixed, data, random = NULL,
                    residual = ~ ar1(row):ar1(col),
                    control = ofe_control()) {
  cl <- match.call()
  if (!inherits(fixed, "formula") || length(fixed) != 3L) {
    stop("`fixed` must be a two-sided formula, e.g. yield ~ treat.",
         call. = FALSE)
  }
  data <- as.data.frame(data)

  vars <- unique(c(all.vars(fixed), all.vars(residual),
                   if (!is.null(random)) all.vars(random)))
  missing_cols <- setdiff(vars, names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  keep <- stats::complete.cases(data[, vars, drop = FALSE])
  n_drop <- sum(!keep)
  dat <- droplevels(data[keep, , drop = FALSE])
  if (nrow(dat) == 0L) stop("No complete observations.", call. = FALSE)
  if (nrow(dat) > control$max_n) {
    stop("This engine factorises an ", nrow(dat), " x ", nrow(dat),
         " matrix at every iteration, above the `max_n` limit of ",
         control$max_n, ". Aggregate to a coarser lattice with ",
         "grid_dense_layer(), or raise ofe_control(max_n = ) knowingly.",
         call. = FALSE)
  }

  mf <- stats::model.frame(fixed, data = dat)
  y <- stats::model.response(mf)
  # A factor that has collapsed to one level (a zone partition that found no
  # zones, a treatment with one rate left after filtering) makes model.matrix
  # fail deep inside contrasts<-; say what actually happened instead.
  flat <- vapply(mf, function(v) is.factor(v) && nlevels(droplevels(v)) < 2L,
                 logical(1))
  if (any(flat)) {
    stop("Fixed-effect factor(s) with a single level: ",
         paste(names(mf)[flat], collapse = ", "),
         ". Drop them from `fixed`.", call. = FALSE)
  }
  mt <- stats::terms(mf)
  X_full <- stats::model.matrix(mt, mf)
  qrX <- qr(X_full)
  aliased <- character(0)
  if (qrX$rank < ncol(X_full)) {
    keep_cols <- sort(qrX$pivot[seq_len(qrX$rank)])
    aliased <- colnames(X_full)[-keep_cols]
    X <- X_full[, keep_cols, drop = FALSE]
    assign_vec <- attr(X_full, "assign")[keep_cols]
  } else {
    X <- X_full
    assign_vec <- attr(X_full, "assign")
  }
  n <- nrow(X)
  p <- ncol(X)
  if (n <= p) {
    stop("Only ", n, " observations for ", p, " fixed-effect parameters.",
         call. = FALSE)
  }

  spec <- .build_residual_spec(residual, dat)
  ptab <- spec$ptab
  scales <- rep(1, nrow(ptab))
  for (sec in spec$sections) {
    for (cp in sec$comps) {
      if (identical(cp$fn, "exp")) scales[cp$par] <- cp$scale
    }
  }

  ZZt <- list()
  ran_par <- integer(0)
  ran_labels <- character(0)
  if (!is.null(random)) {
    if (!inherits(random, "formula") || length(random) != 2L) {
      stop("`random` must be a one-sided formula, e.g. ~ rep.", call. = FALSE)
    }
    ran_labels <- attr(stats::terms(random), "term.labels")
    if (length(ran_labels) == 0L) {
      stop("`random` names no terms.", call. = FALSE)
    }
    fdat <- dat
    for (v in all.vars(random)) fdat[[v]] <- factor(fdat[[v]])
    for (k in seq_along(ran_labels)) {
      Zk <- stats::model.matrix(
        stats::as.formula(paste("~ 0 +", ran_labels[k])), data = fdat)
      ZZt[[k]] <- tcrossprod(Zk)
      ptab <- rbind(ptab, data.frame(name = ran_labels[k], type = "ratio",
                                     stringsAsFactors = FALSE))
      scales <- c(scales, 1)
      ran_par <- c(ran_par, nrow(ptab))
    }
  }

  cache <- list(y = y, X = X, n = n, p = p, spec = spec,
                n_random = length(ZZt), ZZt = ZZt, ran_par = ran_par)

  start_nat <- vapply(seq_len(nrow(ptab)), function(i) {
    switch(ptab$type[i], cor = 0.5, range = scales[i], var = 1, ratio = 0.2)
  }, numeric(1))
  names(start_nat) <- ptab$name
  if (!is.null(control$start)) {
    hit <- intersect(names(control$start), names(start_nat))
    start_nat[hit] <- control$start[hit]
  }
  t0 <- .ofe_to_raw(start_nat, ptab, scales)

  obj <- function(t) {
    v <- .ofe_eval(.ofe_to_nat(t, ptab, scales), cache)
    if (control$trace) message(sprintf("  REML logLik = %12.4f", -v))
    v
  }

  converged <- TRUE
  if (nrow(ptab) == 0L) {
    t_hat <- numeric(0)
  } else {
    # Nelder-Mead first to get out of a bad start, then BFGS to polish. With a
    # single parameter Nelder-Mead is documented as unreliable, so skip it.
    t_start <- if (length(t0) > 1L) {
      stats::optim(t0, obj, method = "Nelder-Mead",
                   control = list(maxit = control$maxit,
                                  reltol = control$tol))$par
    } else {
      t0
    }
    o2 <- stats::optim(t_start, obj, method = "BFGS",
                       control = list(maxit = control$maxit,
                                      reltol = control$tol))
    t_hat <- o2$par
    converged <- o2$convergence == 0L
  }

  theta <- .ofe_to_nat(t_hat, ptab, scales)
  names(theta) <- ptab$name
  res <- .ofe_eval(theta, cache, full = TRUE)
  if (is.null(res)) {
    stop("The covariance matrix is not positive definite at the optimum. ",
         "Try different starting values via ofe_control(start = ).",
         call. = FALSE)
  }

  se_theta <- rep(NA_real_, length(theta))
  if (control$se && length(theta) > 0L) {
    hess <- tryCatch(stats::optimHess(t_hat, obj), error = function(e) NULL)
    if (!is.null(hess)) {
      vc <- tryCatch(solve(hess), error = function(e) NULL)
      if (!is.null(vc) && all(diag(vc) > 0)) {
        jac <- .ofe_jacobian(t_hat, ptab, scales)
        se_theta <- sqrt(diag(vc)) * abs(jac)
      }
    }
  }

  beta <- res$beta
  names(beta) <- colnames(X)
  dimnames(res$vcov) <- list(colnames(X), colnames(X))
  fitted <- drop(X %*% beta)

  varcomp <- data.frame(
    component = c(ptab$name, "sigma2"),
    type = c(ptab$type, "var"),
    estimate = c(unname(theta), res$sigma2),
    std.error = c(se_theta, res$sigma2 * sqrt(2 / (n - p))),
    stringsAsFactors = FALSE
  )
  # Variance ratios are estimated relative to sigma^2; show the variance too,
  # which is what a reader compares against the residual.
  varcomp$variance <- ifelse(varcomp$type %in% c("var", "ratio"),
                             varcomp$estimate * res$sigma2, NA_real_)
  varcomp$variance[varcomp$component == "sigma2"] <- res$sigma2

  structure(list(
    call = cl, fixed = fixed, random = random, residual = residual,
    coefficients = beta, vcov = res$vcov, sigma2 = res$sigma2,
    theta = theta, varcomp = varcomp, loglik = res$loglik,
    fitted = fitted, residuals = as.numeric(y) - fitted,
    n = n, p = p, n_dropped = n_drop, converged = converged,
    aliased = aliased, assign = assign_vec, terms = mt,
    xlevels = stats::.getXlevels(mt, mf), data = dat,
    sections = vapply(spec$sections, function(s) s$tag, character(1))
  ), class = "ofe_fit")
}
