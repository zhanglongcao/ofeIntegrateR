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
#' @param cor_structure Character; `"independent"` (default) fits
#'   `diag(layer):id(unit)`, i.e. separate spatial variances for the dense
#'   and point-source layers with no cross-covariance — robust when point
#'   samples are sparse. `"unstructured"` fits `us(layer):id(unit)`, which
#'   additionally estimates the spatial cross-covariance between layers
#'   (how strongly point-source spatial pattern informs the dense
#'   response), but is often unidentifiable below `min_point_n_for_us`
#'   observations.
#' @param min_point_n_for_us Integer; minimum number of non-missing
#'   point-source observations below which a warning is issued when
#'   `cor_structure = "unstructured"`.
#' @param maxit Maximum number of asreml iterations.
#' @param ai_sing Logical; if `TRUE` (default), temporarily set the
#'   session-level `asreml.options(ai.sing = TRUE)` for the duration of
#'   this fit (restored to its previous value on exit), so that
#'   singularities in the Average Information matrix do not abort
#'   estimation. In testing, the layer-sharing random effect (even with
#'   `cor_structure = "independent"`) can trigger this with sparse
#'   point-source data, so it is enabled by default.
#' @param ... Additional arguments passed to [asreml::asreml()].
#'
#' @return A fitted `asreml` model object. Treatment contrasts for the
#'   dense layer can be retrieved with [extract_fixed_effects()] (look for
#'   terms matching `at(layer, "dense"):<treat>`).
#'
#' @details
#' Residuals are modelled with `dsum(~units | layer)`, i.e. separate i.i.d.
#' residual variances per layer. An `ar1(row):ar1(col)` residual was found
#' to be unidentifiable for the point-source layer in testing, because a
#' handful of sparse, irregularly placed observations cannot support a
#' full spatial autocorrelation structure on their own; the shared random
#' effect (`diag(layer)` or `us(layer)`) is what carries spatial
#' information between layers instead.
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
#' @export
fit_integrated_joint <- function(grid, point_samples,
                                  response_dense, response_point,
                                  treat = "treat", row = "row", col = "col",
                                  cor_structure = c("independent", "unstructured"),
                                  min_point_n_for_us = 30,
                                  maxit = 60, ai_sing = TRUE, ...) {
  cor_structure <- match.arg(cor_structure)
  check_asreml()

  n_point <- sum(!is.na(point_samples[[response_point]]))
  if (cor_structure == "unstructured" && n_point < min_point_n_for_us) {
    warning(
      "Only ", n_point, " point-source observations, fewer than ",
      "`min_point_n_for_us` = ", min_point_n_for_us, ". The unstructured ",
      "cross-covariance (us(layer)) is often unidentifiable with sparse ",
      "point data; consider cor_structure = \"independent\".", call. = FALSE
    )
  }

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
