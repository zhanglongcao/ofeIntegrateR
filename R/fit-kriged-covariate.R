#' Fit a treatment model using a kriged point-source covariate
#'
#' Implements the "kriged-covariate" data integration strategy: a sparse
#' point-source variable (already kriged onto the trial grid, e.g. via
#' [krige_point_samples()]) is included as a fixed covariate alongside the
#' treatment factor, with an AR1xAR1 spatial residual structure (when
#' `engine = "asreml"`) to absorb any remaining spatial autocorrelation in
#' the dense response.
#'
#' @param data Data frame containing the response, treatment factor, kriged
#'   covariate, and row/column position columns (one row per grid cell).
#' @param response Character; name of the dense response column (e.g.
#'   yield).
#' @param treat Character; name of the treatment factor column.
#' @param covariate Character; name of the kriged point-source covariate
#'   column (e.g. the `*_kriged` column produced by [krige_point_samples()]).
#' @param row,col Character; names of the row/column position columns used
#'   to build the spatial residual structure. Ignored when `engine = "lm"`.
#' @param engine Character; `"asreml"` (default) fits
#'   `response ~ treat + covariate` with an `ar1(row):ar1(col)` residual via
#'   [asreml::asreml()] — requires a licensed copy of asreml-R. `"gls"` fits
#'   the same fixed-effects model via [nlme::gls()] with an exponential
#'   spatial correlation structure ([nlme::corExp()] on `row`/`col`) — an
#'   open-source (CRAN-only, no licence required) alternative that still
#'   accounts for residual spatial autocorrelation. `"lm"` fits the same
#'   fixed-effects model via [stats::lm()] with no spatial residual
#'   structure at all — the simplest fallback, but does not account for
#'   residual spatial autocorrelation.
#' @param ... Additional arguments passed to [asreml::asreml()] (e.g.
#'   `maxit`), [nlme::gls()], or [stats::lm()].
#'
#' @return The fitted model object (class `asreml`, `gls`, or `lm`). Use
#'   [extract_fixed_effects()] to retrieve a tidy table of fixed-effect
#'   estimates, including the treatment contrasts.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 12, seed = 1)
#' krieged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
#' fit_lm <- fit_integrated_kriged(krieged, response = "dense_response",
#'                                  treat = "treat", covariate = "point_obs_kriged",
#'                                  row = "row", col = "col", engine = "lm")
#' extract_fixed_effects(fit_lm)
#'
#' if (requireNamespace("asreml", quietly = TRUE)) {
#'   fit_asr <- fit_integrated_kriged(krieged, response = "dense_response",
#'                                     treat = "treat", covariate = "point_obs_kriged",
#'                                     row = "row", col = "col", engine = "asreml")
#'   extract_fixed_effects(fit_asr)
#' }
#'
#' @export
fit_integrated_kriged <- function(data, response, treat, covariate,
                                   row = "row", col = "col",
                                   engine = c("asreml", "gls", "lm"), ...) {
  engine <- match.arg(engine)
  needed <- c(response, treat, covariate)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  fixed <- stats::as.formula(paste(response, "~", treat, "+", covariate))

  if (engine == "lm") {
    return(stats::lm(fixed, data = data, ...))
  }

  if (engine == "gls") {
    data$.row_num <- as.numeric(as.character(data[[row]]))
    data$.col_num <- as.numeric(as.character(data[[col]]))
    return(nlme::gls(
      fixed,
      data = data,
      correlation = nlme::corExp(form = ~ .row_num + .col_num, nugget = TRUE),
      ...
    ))
  }

  check_asreml()
  data[[row]] <- as.factor(data[[row]])
  data[[col]] <- as.factor(data[[col]])
  residual <- stats::as.formula(paste0("~ ar1(", row, "):ar1(", col, ")"))

  asreml::asreml(
    fixed = fixed,
    residual = residual,
    data = data,
    trace = FALSE,
    ...
  )
}

#' Extract a tidy table of fixed-effect estimates
#'
#' Works for `lm`, `gls`, `asreml`, and `mmer` model objects, so the same
#' downstream code can summarise treatment contrasts regardless of which
#' `engine` was used in [fit_integrated_kriged()] or [fit_integrated_joint()].
#'
#' @param model A fitted model of class `lm`, `gls`, `asreml`, or `mmer`.
#'
#' @return A data frame with columns `term`, `estimate`, `se`.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 12, seed = 1)
#' krieged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
#' fit_lm <- fit_integrated_kriged(krieged, response = "dense_response",
#'                                  treat = "treat", covariate = "point_obs_kriged",
#'                                  row = "row", col = "col", engine = "lm")
#' extract_fixed_effects(fit_lm)
#'
#' @export
extract_fixed_effects <- function(model) {
  if (inherits(model, "asreml")) {
    cf <- summary(model, coef = TRUE)$coef.fixed
    return(data.frame(
      term = rownames(cf),
      estimate = cf[, "solution"],
      se = cf[, "std error"],
      row.names = NULL
    ))
  }
  if (inherits(model, "gls")) {
    cf <- summary(model)$tTable
    return(data.frame(
      term = rownames(cf),
      estimate = cf[, "Value"],
      se = cf[, "Std.Error"],
      row.names = NULL
    ))
  }
  if (inherits(model, "mmer")) {
    cf <- model$Beta
    return(data.frame(
      term = paste(cf$Trait, cf$Effect, sep = ":"),
      estimate = cf$Estimate,
      se = sqrt(diag(model$VarBeta)),
      row.names = NULL
    ))
  }
  if (inherits(model, "lm")) {
    cf <- summary(model)$coefficients
    return(data.frame(
      term = rownames(cf),
      estimate = cf[, "Estimate"],
      se = cf[, "Std. Error"],
      row.names = NULL
    ))
  }
  stop("Unsupported model class: ", paste(class(model), collapse = "/"),
       ". Expected 'lm', 'gls', 'asreml', or 'mmer'.", call. = FALSE)
}
