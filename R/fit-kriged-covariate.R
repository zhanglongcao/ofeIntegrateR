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
#' @param covariate Character vector; name(s) of the kriged point-source
#'   covariate column(s) (e.g. the `*_kriged` columns produced by
#'   [krige_point_samples()]). Supply several names to integrate more than one
#'   point variable at once, such as soil nitrogen and phosphorus. Use
#'   `NULL` to fit the treatment model with no point-source covariate at all:
#'   this is the baseline that an integrated analysis has to beat, since it
#'   uses only the dense layer and costs no sampling.
#' @param row,col Character; names of the row/column position columns used
#'   to build the spatial residual structure. Ignored when `engine = "lm"`.
#' @param engine Character; `"asreml"` (default) fits
#'   `response ~ treat + covariate` with an `ar1(row):ar1(col)` residual via
#'   [asreml::asreml()] — requires a licensed copy of asreml-R. `"lme"` fits
#'   random effects **and** an exponential spatial correlation via
#'   [nlme::lme()], which is the open-source counterpart to the asreml fit for
#'   a replicated trial. `"gls"` fits
#'   the same fixed-effects model via [nlme::gls()] with an exponential
#'   spatial correlation structure ([nlme::corExp()] on `row`/`col`) — an
#'   open-source (CRAN-only, no licence required) alternative that still
#'   accounts for residual spatial autocorrelation. `"lm"` fits the same
#'   fixed-effects model via [stats::lm()] with no spatial residual
#'   structure at all — the simplest fallback, but does not account for
#'   residual spatial autocorrelation.
#' @param random Optional one-sided formula of random effects, such as
#'   `~ rep` or `~ block`. Real strip trials are replicated, so this is usually
#'   needed. Supported by `engine = "asreml"` and `engine = "lme"`; the `"gls"`
#'   and `"lm"` engines cannot fit random effects and raise an error rather
#'   than ignoring the argument.
#' @param ... Additional arguments passed to [asreml::asreml()] (e.g.
#'   `maxit`), [nlme::lme()], [nlme::gls()], or [stats::lm()].
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
# asreml takes `random = ~ rep`; nlme::lme wants a grouping formula,
# `~ 1 | rep`. Accepting the asreml spelling on both engines is the point of a
# common interface, so translate rather than making the caller remember which
# is which. A formula that already names a grouping factor is passed through.
as_lme_random <- function(random) {
  if (inherits(random, "formula") && length(random) == 2L) {
    rhs <- random[[2]]
    if (!(is.call(rhs) && identical(as.character(rhs[[1]]), "|"))) {
      terms_chr <- attr(stats::terms(random), "term.labels")
      if (length(terms_chr) == 0L) {
        stop("`random` names no grouping factor.", call. = FALSE)
      }
      if (length(terms_chr) > 1L) {
        stop("`engine = \"lme\"` needs the nesting made explicit for more ",
             "than one grouping factor, e.g. random = ~ 1 | block/plot ",
             "instead of ~ block + plot.", call. = FALSE)
      }
      return(stats::as.formula(paste("~ 1 |", terms_chr)))
    }
  }
  random
}

fit_integrated_kriged <- function(data, response, treat, covariate = NULL,
                                   random = NULL,
                                   row = "row", col = "col",
                                   engine = c("asreml", "lme", "gls", "lm"),
                                   ...) {
  engine <- match.arg(engine)
  if (!is.null(random) && engine %in% c("gls", "lm")) {
    stop("`engine = \"", engine, "\"` cannot fit random effects. Use ",
         "engine = \"lme\" for an open-source fit with both random effects ",
         "and a spatial correlation structure, or engine = \"asreml\" if ",
         "licensed.", call. = FALSE)
  }
  needed <- c(response, treat, covariate)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  # `covariate` may name several kriged surfaces (e.g. soil N and P). Collapse
  # them into one right-hand side: paste() vectorises, so without `collapse` the
  # formula would be built from the first covariate only and the rest silently
  # dropped. With no covariate this is the dense-layer-only baseline.
  rhs <- if (length(covariate) == 0L) {
    treat
  } else {
    paste(treat, "+", paste(covariate, collapse = " + "))
  }
  fixed <- stats::as.formula(paste(response, "~", rhs))

  if (engine == "lm") {
    return(stats::lm(fixed, data = data, ...))
  }

  if (engine == "lme") {
    # lme carries both a random-effects structure and a spatial correlation on
    # the residual, so it is the open-source counterpart to the asreml fit for
    # a replicated trial. gls has the correlation but no random effects.
    if (is.null(random)) {
      stop("`engine = \"lme\"` needs a `random` formula, e.g. random = ~ rep. ",
           "With no random effects use engine = \"gls\", which fits the same ",
           "spatial correlation structure.", call. = FALSE)
    }
    data$.row_num <- as.numeric(as.character(data[[row]]))
    data$.col_num <- as.numeric(as.character(data[[col]]))
    dots <- list(...)
    if (is.null(dots$na.action)) dots$na.action <- stats::na.omit
    return(do.call(nlme::lme, c(
      list(fixed = fixed,
           data = data,
           random = as_lme_random(random),
           correlation = nlme::corExp(form = ~ .row_num + .col_num,
                                      nugget = TRUE)),
      dots
    )))
  }

  if (engine == "gls") {
    data$.row_num <- as.numeric(as.character(data[[row]]))
    data$.col_num <- as.numeric(as.character(data[[col]]))
    # Gridded trial data routinely contain empty cells and cells excluded for
    # straddling a treatment boundary. `gls()` stops on missing values rather
    # than dropping them, so omit incomplete rows unless the caller has said
    # how to handle them. `gls()` uses actual coordinates rather than a
    # complete lattice, so dropping rows is safe here (it is not for the
    # separable AR1 residual used by the asreml engine, which keeps NA rows).
    dots <- list(...)
    if (is.null(dots$na.action)) dots$na.action <- stats::na.omit
    return(do.call(nlme::gls, c(
      list(model = fixed,
           data = data,
           correlation = nlme::corExp(form = ~ .row_num + .col_num,
                                      nugget = TRUE)),
      dots
    )))
  }

  check_asreml()
  data[[row]] <- as.factor(data[[row]])
  data[[col]] <- as.factor(data[[col]])
  residual <- stats::as.formula(paste0("~ ar1(", row, "):ar1(", col, ")"))

  asreml_args <- list(fixed = fixed, residual = residual, data = data,
                      trace = FALSE)
  if (!is.null(random)) asreml_args$random <- random
  do.call(asreml::asreml, c(asreml_args, list(...)))
}

#' Extract a tidy table of fixed-effect estimates
#'
#' Works for `lm`, `gls`, `lme`, `asreml`, and `mmer` model objects, so the same
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
  if (inherits(model, "lme")) {
    cf <- summary(model)$tTable
    return(data.frame(
      term = rownames(cf),
      estimate = cf[, "Value"],
      se = cf[, "Std.Error"],
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
