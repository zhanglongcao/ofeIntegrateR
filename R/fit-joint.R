#' Build a long-format combined data frame for the joint model
#'
#' Stacks the dense response and the (sparse) point-source response into
#' two "layers" sharing the same grid locations, as required by
#' [fit_integrated_joint()]. The point-source layer carries `NA` at every
#' location not in `point_samples`.
#'
#' @inheritParams fit_integrated_joint
#' @keywords internal
#' @noRd
build_long_layers <- function(grid, point_samples, response_dense, response_point,
                               treat, row, col) {
  dense_layer <- grid
  dense_layer$response <- grid[[response_dense]]
  dense_layer$layer <- "dense"
  dense_layer$unit <- paste(grid[[row]], grid[[col]], sep = "_")

  point_layer <- grid
  point_layer$response <- NA_real_
  point_layer$layer <- "point"
  point_layer$unit <- paste(grid[[row]], grid[[col]], sep = "_")

  for (i in seq_len(nrow(point_samples))) {
    hit <- which(point_layer[[row]] == point_samples[[row]][i] &
                   point_layer[[col]] == point_samples[[col]][i])
    if (length(hit)) {
      point_layer$response[hit] <- point_samples[[response_point]][i]
    }
  }

  keep <- c("response", "layer", "unit", row, col, treat)
  combined <- rbind(dense_layer[, keep], point_layer[, keep])
  combined$layer <- factor(combined$layer, levels = c("dense", "point"))
  combined[[row]] <- factor(combined[[row]])
  combined[[col]] <- factor(combined[[col]])
  combined$unit <- factor(combined$unit)
  combined
}

#' Fit the joint bivariate model with sommer (open-source fallback)
#'
#' @inheritParams fit_integrated_joint
#' @keywords internal
#' @noRd
fit_joint_sommer <- function(grid, point_samples, response_dense, response_point,
                              treat, row, col, cor_structure, ...) {
  check_sommer()
  wide <- grid
  wide$unit <- factor(paste(wide[[row]], wide[[col]], sep = "_"))
  wide[[response_point]] <- NA_real_
  for (i in seq_len(nrow(point_samples))) {
    hit <- which(wide[[row]] == point_samples[[row]][i] &
                   wide[[col]] == point_samples[[col]][i])
    if (length(hit)) {
      wide[[response_point]][hit] <- point_samples[[response_point]][i]
    }
  }

  fixed <- stats::as.formula(
    paste0("cbind(`", response_dense, "`, `", response_point, "`) ~ ", treat)
  )
  random <- if (cor_structure == "independent") {
    ~ sommer::vsr(unit, Gtc = base::diag(2))
  } else {
    ~ sommer::vsr(unit, Gtc = sommer::unsm(2))
  }
  rcov <- ~ sommer::vsr(units, Gtc = base::diag(2))

  sommer::mmer(
    fixed = fixed,
    random = random,
    rcov = rcov,
    data = wide,
    naMethodY = "include",
    verbose = FALSE,
    ...
  )
}

#' Fit a joint bivariate spatial model for dense and point-source data
#'
#' Implements the "joint model" data integration strategy: rather than
#' pre-kriging the point-source variable and using it as a fixed covariate
#' (see [fit_integrated_kriged()]), the dense response and the point-source
#' measurements are fitted **jointly** as two stacked response "layers"
#' sharing a spatial random effect at each grid location. This lets sparse
#' point samples directly inform the dense response's spatial structure
#' (and vice versa) without a separate kriging step.
#'
#' @param grid Data frame, one row per grid location, with the dense
#'   response, treatment factor, and row/column position columns (as
#'   produced by [simulate_ofe_trial()]'s `grid` element).
#' @param point_samples Data frame of sparse point-source observations with
#'   row/column position columns and the point-source value column (as
#'   produced by [simulate_ofe_trial()]'s `point_samples` element).
#' @param response_dense Character; name of the dense response column in
#'   `grid` (e.g. yield).
#' @param response_point Character; name of the point-source value column
#'   in `point_samples` (e.g. soil, tissue, or disease measurement).
#' @param treat Character; name of the treatment factor column in `grid`.
#'   The treatment effect is estimated from the dense layer only.
#' @param row,col Character; names of the row/column position columns,
#'   shared between `grid` and `point_samples`.
#' @param engine Character; `"asreml"` (default) fits the model via
#'   `asreml::asreml()` — requires a licensed copy of asreml-R. `"sommer"`
#'   fits an open-source (CRAN-only, no licence required) analogue via
#'   `sommer::mmer()`: the dense and point-source responses are modelled as
#'   two traits in a multi-trait model sharing a random effect per grid
#'   location (`vsr(unit, Gtc = ...)`), with the treatment effect estimated
#'   separately for each trait (only the dense-trait estimate is
#'   meaningful; see Details).
#' @param cor_structure Character; `"independent"` (default) fits
#'   `diag(layer):id(unit)` (asreml) or `vsr(unit, Gtc = diag(2))` (sommer),
#'   i.e. separate spatial variances for the dense and point-source layers
#'   with no cross-covariance — robust when point samples are sparse.
#'   `"unstructured"` fits `us(layer):id(unit)` (asreml) or
#'   `vsr(unit, Gtc = unsm(2))` (sommer), which additionally estimates the
#'   spatial cross-covariance between layers (how strongly point-source
#'   spatial pattern informs the dense response), but is often
#'   unidentifiable below `min_point_n_for_us` observations.
#' @param min_point_n_for_us Integer; minimum number of non-missing
#'   point-source observations below which a warning is issued when
#'   `cor_structure = "unstructured"`.
#' @param maxit Maximum number of asreml iterations. Ignored when
#'   `engine = "sommer"`.
#' @param ai_sing Logical; if `TRUE` (default), temporarily set the
#'   session-level `asreml.options(ai.sing = TRUE)` for the duration of
#'   this fit (restored to its previous value on exit), so that
#'   singularities in the Average Information matrix do not abort
#'   estimation. In testing, the layer-sharing random effect (even with
#'   `cor_structure = "independent"`) can trigger this with sparse
#'   point-source data, so it is enabled by default. Ignored when
#'   `engine = "sommer"`.
#' @param ... Additional arguments passed to `asreml::asreml()` or
#'   `sommer::mmer()`.
#'
#' @return A fitted `asreml` or `mmer` model object. Treatment contrasts
#'   for the dense layer can be retrieved with [extract_fixed_effects()]
#'   (look for terms matching `at(layer, "dense"):<treat>` for asreml, or
#'   `<response_dense>:<treat><level>` for sommer).
#'
#' @details
#' For `engine = "asreml"`, residuals are modelled with
#' `dsum(~units | layer)`, i.e. separate i.i.d. residual variances per
#' layer. An `ar1(row):ar1(col)` residual was found to be unidentifiable
#' for the point-source layer in testing, because a handful of sparse,
#' irregularly placed observations cannot support a full spatial
#' autocorrelation structure on their own; the shared random effect
#' (`diag(layer)` or `us(layer)`) is what carries spatial information
#' between layers instead.
#'
#' For `engine = "sommer"`, the same idea is expressed as a multi-trait
#' model (`cbind(dense, point) ~ treat`) with a shared random effect per
#' grid location (`unit`) and heterogeneous, uncorrelated residual
#' variances per trait (`vsr(units, Gtc = diag(2))`). sommer estimates the
#' treatment effect separately for *each* trait by default; the
#' point-source trait's treatment estimate is not meaningful (point
#' samples carry no treatment information of their own) and should be
#' ignored — only the dense-trait estimate is the quantity of interest.
#' Rows are kept even when only one of the two responses is observed via
#' `naMethodY = "include"`. In testing against simulated data, this
#' sommer-based fit reproduced the asreml `cor_structure = "independent"`
#' treatment-contrast estimates almost exactly.
#'
#' @examples
#' if (requireNamespace("asreml", quietly = TRUE)) {
#'   sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 12, seed = 1)
#'   fit <- fit_integrated_joint(sim$grid, sim$point_samples,
#'                                response_dense = "dense_response",
#'                                response_point = "point_obs")
#'   extract_fixed_effects(fit)
#' }
#'
#' if (requireNamespace("sommer", quietly = TRUE)) {
#'   # sommer's multi-trait model needs more point-source observations than
#'   # asreml's to avoid a singular system; a sparser grid like the asreml
#'   # example above is not guaranteed to converge.
#'   sim <- simulate_ofe_trial(n_row = 30, n_col = 21, n_point_samples = 25, seed = 11)
#'   fit <- fit_integrated_joint(sim$grid, sim$point_samples,
#'                                response_dense = "dense_response",
#'                                response_point = "point_obs",
#'                                engine = "sommer")
#'   extract_fixed_effects(fit)
#' }
#'
#' @export
fit_integrated_joint <- function(grid, point_samples,
                                  response_dense, response_point,
                                  treat = "treat", row = "row", col = "col",
                                  engine = c("asreml", "sommer"),
                                  cor_structure = c("independent", "unstructured"),
                                  min_point_n_for_us = 30,
                                  maxit = 60, ai_sing = TRUE, ...) {
  engine <- match.arg(engine)
  cor_structure <- match.arg(cor_structure)

  n_point <- sum(!is.na(point_samples[[response_point]]))
  if (cor_structure == "unstructured" && n_point < min_point_n_for_us) {
    warning(
      "Only ", n_point, " point-source observations, fewer than ",
      "`min_point_n_for_us` = ", min_point_n_for_us, ". The unstructured ",
      "cross-covariance is often unidentifiable with sparse point data; ",
      "consider cor_structure = \"independent\".", call. = FALSE
    )
  }

  if (engine == "sommer") {
    return(fit_joint_sommer(grid, point_samples, response_dense, response_point,
                             treat, row, col, cor_structure, ...))
  }

  check_asreml()
  combined <- build_long_layers(grid, point_samples, response_dense, response_point,
                                 treat, row, col)

  fixed <- stats::as.formula(
    paste0('response ~ layer + at(layer, "dense"):', treat)
  )

  random <- if (cor_structure == "independent") {
    ~ diag(layer):id(unit)
  } else {
    ~ us(layer, init = c(0.8, 0, 0.8)):id(unit)
  }

  if (ai_sing) {
    # ai.sing is a session-level option in asreml-R (>= 4.2), not an
    # argument to asreml(); passing it directly to asreml() is silently
    # ignored. Set it for the duration of this call only.
    old_ai_sing <- asreml::asreml.options()$ai.sing
    on.exit(asreml::asreml.options(ai.sing = old_ai_sing), add = TRUE)
    asreml::asreml.options(ai.sing = TRUE)
  }

  asreml::asreml(
    fixed = fixed,
    random = random,
    residual = ~ dsum(~units | layer),
    data = combined,
    na.action = asreml::na.method(y = "include", x = "include"),
    maxit = maxit,
    trace = FALSE,
    ...
  )
}
