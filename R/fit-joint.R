#' Put the sparse layer onto the dense layer's grid, one column per trait
#'
#' Matches each point sample to its grid cell by `row` and `col`. Several
#' samples in one cell are averaged; samples that fall on no cell are dropped
#' with a warning, since silently losing them would change what the
#' cross-covariance is estimated from.
#'
#' @inheritParams fit_integrated_joint
#' @keywords internal
#' @noRd
build_joint_wide <- function(grid, point_samples, response_dense, response_point,
                             row, col) {
  key_g <- paste(grid[[row]], grid[[col]], sep = "_")
  key_p <- paste(point_samples[[row]], point_samples[[col]], sep = "_")
  val <- point_samples[[response_point]]
  keep <- !is.na(val)
  key_p <- key_p[keep]; val <- val[keep]

  off <- !key_p %in% key_g
  if (any(off)) {
    warning(sum(off), " point sample(s) fall on no cell of `grid` and were ",
            "dropped. Check that `row` and `col` index the same lattice in ",
            "both data frames.", call. = FALSE)
    key_p <- key_p[!off]; val <- val[!off]
  }
  if (anyDuplicated(key_p)) {
    message("Several point samples share a grid cell; their mean is used.")
  }
  cell_mean <- tapply(val, key_p, mean)

  wide <- grid
  wide[[response_point]] <- unname(cell_mean[key_g])
  wide
}

#' Fit the joint bivariate model with sommer (open-source fallback)
#'
#' @inheritParams fit_integrated_joint
#' @keywords internal
#' @noRd
fit_joint_sommer <- function(wide, response_dense, response_point, treat, block,
                             ...) {
  check_sommer()
  # sommer::mmer() cannot leave a response missing: naMethodY = "include" and
  # "include2" impute it with the trait median, and "exclude" drops the cell.
  # For a sparse layer, imputation swamps the few real observations (the
  # layer's treatment effects shrink towards zero and its SEs collapse), and
  # dropping cells discards the dense layer the model exists to use. Either
  # way the answer is wrong, so only a lattice complete in both layers is
  # accepted here.
  n_miss <- sum(!stats::complete.cases(wide[, c(response_dense, response_point)]))
  if (n_miss > 0) {
    stop("engine = \"sommer\" needs both layers observed in every cell, but ",
         n_miss, " cell(s) lack one of them. sommer::mmer() imputes a missing ",
         "response with its median or drops the cell, and either gives wrong ",
         "estimates for a sparse layer. Use engine = \"asreml\".",
         call. = FALSE)
  }
  rhs <- paste(c(treat, block), collapse = " + ")
  fixed <- stats::as.formula(
    paste0("cbind(`", response_dense, "`, `", response_point, "`) ~ ", rhs)
  )
  # units:us(trait) -- the cross-covariance lives in the residual, which is
  # where it is identifiable; sommer has no ar1 x ar1 x us(trait) residual
  sommer::mmer(
    fixed = fixed,
    rcov = ~ sommer::vsr(units, Gtc = sommer::unsm(2)),
    data = wide,
    naMethodY = "include",
    verbose = FALSE,
    ...
  )
}

#' Fit a joint bivariate spatial model for a dense and a sparse layer
#'
#' Implements the "joint model" data integration strategy. The dense layer
#' (e.g. yield-monitor yield, NDVI) and the sparse layer (e.g. soil cores,
#' hand-harvest cuts, grain protein) are modelled as two correlated traits,
#' each with its own treatment effects, on the same lattice. The sparse layer
#' is missing at every cell that was not sampled.
#'
#' Precision is gained through the conditional distribution. If the two
#' layers' residuals at a cell are bivariate normal with correlation
#' \eqn{\rho}, the residual of one given the other has variance
#' \eqn{\sigma^2(1-\rho^2)}. The model therefore borrows strength only through
#' the residual cross-covariance, and the borrowing mostly flows *from* the
#' dense layer *to* the sparse one. The practical consequence:
#'
#' * **Sparse response, dense auxiliary layer** -- for example protein or
#'   hand cuts with NDVI or the yield monitor alongside. This is where the model
#'   earns its keep: in simulation, the treatment-contrast SD of the sparse
#'   response fell by 16--38% relative to analysing it alone, with nominal
#'   coverage.
#' * **Dense response, sparse covariate** -- for example yield with 30 soil
#'   cores. The dense layer already carries almost all the information, and
#'   the gain over the free spatial baseline was within a few percent. Here
#'   [fit_integrated_kriged()] does as well and is cheaper.
#'
#' Both layers are given treatment effects. That is required when the
#' auxiliary layer responds to treatment (NDVI after a nitrogen treatment,
#' say): conditioning is then on the deviations from each layer's own
#' treatment means, so the auxiliary layer cannot absorb the treatment effect
#' of the response. This is the model-based counterpart of centring a
#' covariate within plots (Piepho et al.), and it needs no centring by hand.
#' For a pre-treatment layer the extra treatment terms cost a few degrees of
#' freedom and nothing else.
#'
#' @param grid Data frame, one row per cell of a **complete** row-by-column
#'   lattice (as returned by [grid_dense_layer()] with `keep_empty = TRUE`),
#'   holding the dense layer, the treatment factor and the position columns.
#' @param point_samples Data frame of the sparse layer, with `row`/`col`
#'   position columns on the same lattice as `grid` and the sparse value
#'   column. Several samples in one cell are averaged.
#' @param response_dense Character; name of the dense layer's column in
#'   `grid`.
#' @param response_point Character; name of the sparse layer's column in
#'   `point_samples`.
#' @param treat Character; name of the treatment factor column in `grid`.
#' @param block Optional character; name of a blocking factor in `grid`
#'   (e.g. replicate), fitted as a fixed effect within each trait.
#' @param row,col Character; names of the row/column position columns,
#'   shared between `grid` and `point_samples`.
#' @param engine `"asreml"` (default; requires a licence) or `"sommer"`
#'   (open source). sommer accepts only a lattice on which both layers are
#'   observed in every cell: `sommer::mmer()` imputes a missing response with
#'   its median or drops the cell, which gives wrong estimates for a sparse
#'   layer, so a sparse layer needs asreml.
#' @param residual `"ar1"` fits the separable bivariate spatial residual
#'   `ar1(row):ar1(col):us(trait)`: the dense layer carries the AR1
#'   parameters, and `us(trait)` carries the cross-covariance. `"id"` fits
#'   `units:us(trait)`, the non-spatial version, which is equivalent to an
#'   analysis of covariance on the within-plot deviations. Defaults to `"ar1"`
#'   for asreml and `"id"` for sommer, which cannot fit the spatial form.
#' @param min_overlap Integer; a warning is issued when fewer cells than this
#'   carry both layers, since those cells are all the cross-covariance is
#'   estimated from.
#' @param maxit Maximum number of asreml iterations. Ignored for sommer.
#' @param ai_sing Logical; if `TRUE` (default), set
#'   `asreml.options(ai.sing = TRUE)` for the duration of the fit (restored on
#'   exit). Ignored for sommer.
#' @param ... Further arguments passed to `asreml::asreml()` or
#'   `sommer::mmer()`.
#'
#' @return A fitted `asreml` or `mmer` object. Read the contrasts with
#'   [extract_fixed_effects()]. With asreml the treatment terms are named
#'   `at(trait, <response>):<treat>_<level>`; with sommer,
#'   `<response>:<treat><level>`. The first treatment level is the reference,
#'   so each estimate is a contrast against it, and its standard error is the
#'   SED.
#'
#' @details
#' Earlier versions fitted a shared `diag(layer):id(unit)` or
#' `us(layer):id(unit)` random effect with independent residuals per layer.
#' With one observation per cell and layer, that unit effect is confounded
#' with the residual. The `diag` form therefore had no cross-covariance at all,
#' and reproduced a univariate i.i.d. analysis of the dense layer exactly,
#' without its spatial structure. The `us` form converged in under 11% of
#' simulated trials. Both are removed.
#'
#' A plot-level cross-covariance (`plot:us(trait)`, as in a split-plot
#' bivariate model) is deliberately not offered. Once treatment and block are
#' fitted, strip trials leave almost no between-strip variation, the plot-level
#' covariance sits on the boundary, and the fit converged in about half of
#' simulated trials without improving precision.
#'
#' @references
#' Piepho, H.-P. et al. (2026). Improving yield estimates in on-farm
#' experiments using vegetation indices as covariates. Manuscript under review.
#'
#' @seealso [fit_integrated_kriged()] for the dense-response, sparse-covariate
#'   case.
#'
#' @examples
#' if (requireNamespace("asreml", quietly = TRUE)) {
#'   sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 30, seed = 1)
#'   fit <- fit_integrated_joint(sim$grid, sim$point_samples,
#'                               response_dense = "dense_response",
#'                               response_point = "point_obs")
#'   fx <- extract_fixed_effects(fit)
#'   fx[grepl("dense_response):treat", fx$term), ]
#' }
#'
#' # sommer fits the non-spatial model when both layers cover every cell,
#' # e.g. a yield map and an NDVI map on the same lattice
#' if (requireNamespace("sommer", quietly = TRUE)) {
#'   sim <- simulate_ofe_trial(n_row = 20, n_col = 12, seed = 11)
#'   both <- sim$grid[, c("row", "col", "point_true")]
#'   fit <- fit_integrated_joint(sim$grid, both,
#'                               response_dense = "dense_response",
#'                               response_point = "point_true",
#'                               engine = "sommer", min_overlap = 0)
#'   extract_fixed_effects(fit)
#' }
#'
#' @export
fit_integrated_joint <- function(grid, point_samples,
                                 response_dense, response_point,
                                 treat = "treat", block = NULL,
                                 row = "row", col = "col",
                                 engine = c("asreml", "sommer"),
                                 residual = NULL,
                                 min_overlap = 10,
                                 maxit = 60, ai_sing = TRUE, ...) {
  engine <- match.arg(engine)
  dots <- list(...)
  if (any(c("cor_structure", "min_point_n_for_us") %in% names(dots))) {
    stop("`cor_structure` and `min_point_n_for_us` have been removed. The ",
         "cross-covariance is now always estimated, in the residual, where ",
         "it is identifiable. See ?fit_integrated_joint, Details.",
         call. = FALSE)
  }
  if (is.null(residual)) residual <- if (engine == "asreml") "ar1" else "id"
  residual <- match.arg(residual, c("ar1", "id"))
  if (engine == "sommer" && residual == "ar1") {
    stop("sommer cannot fit the ar1(row):ar1(col):us(trait) residual. Use ",
         "residual = \"id\", or engine = \"asreml\".", call. = FALSE)
  }

  needed <- c(response_dense, treat, block, row, col)
  miss <- setdiff(needed, names(grid))
  if (length(miss)) {
    stop("`grid` has no column(s): ", paste(miss, collapse = ", "), ".",
         call. = FALSE)
  }
  miss <- setdiff(c(response_point, row, col), names(point_samples))
  if (length(miss)) {
    stop("`point_samples` has no column(s): ", paste(miss, collapse = ", "),
         ".", call. = FALSE)
  }
  resp <- c(response_dense, response_point)
  if (response_dense == response_point) {
    stop("`response_dense` and `response_point` must be different names.",
         call. = FALSE)
  }
  if (!all(make.names(resp) == resp)) {
    stop("Response names must be syntactic R names (they become trait ",
         "levels in the model): ", paste(resp, collapse = ", "), ".",
         call. = FALSE)
  }

  wide <- build_joint_wide(grid, point_samples, response_dense, response_point,
                           row, col)
  wide[[treat]] <- factor(wide[[treat]])
  if (!is.null(block)) wide[[block]] <- factor(wide[[block]])

  n_overlap <- sum(!is.na(wide[[response_dense]]) &
                     !is.na(wide[[response_point]]))
  if (n_overlap < min_overlap) {
    warning("Only ", n_overlap, " cell(s) carry both layers. The ",
            "cross-covariance is estimated from these cells alone and will ",
            "be poorly determined.", call. = FALSE)
  }

  if (engine == "sommer") {
    return(fit_joint_sommer(wide, response_dense, response_point, treat, block,
                            ...))
  }

  check_asreml()
  if (residual == "ar1") {
    cells <- paste(wide[[row]], wide[[col]])
    if (anyDuplicated(cells) ||
        length(cells) != length(unique(wide[[row]])) * length(unique(wide[[col]]))) {
      stop("residual = \"ar1\" needs a complete row-by-column lattice with ",
           "one row per cell. Build `grid` with grid_dense_layer(keep_empty = ",
           "TRUE), or use residual = \"id\".", call. = FALSE)
    }
    wide <- wide[order(wide[[row]], wide[[col]]), ]
  }
  wide[[row]] <- factor(wide[[row]])
  wide[[col]] <- factor(wide[[col]])

  trt_terms <- paste0('at(trait, "', resp, '"):', treat)
  blk_terms <- if (is.null(block)) character(0) else paste0("trait:", block)
  fixed <- stats::as.formula(paste0(
    "cbind(", resp[1], ", ", resp[2], ") ~ ",
    paste(c("trait", blk_terms, trt_terms), collapse = " + ")))
  resid <- if (residual == "ar1") {
    stats::as.formula(paste0("~ ar1(", row, "):ar1(", col, "):us(trait)"))
  } else {
    ~ units:us(trait)
  }

  if (ai_sing) {
    # ai.sing is a session-level option in asreml-R (>= 4.2), not an
    # argument to asreml(); passing it directly to asreml() is silently
    # ignored. Set it for the duration of this call only.
    old_ai_sing <- asreml::asreml.options()$ai.sing
    on.exit(asreml::asreml.options(ai.sing = old_ai_sing), add = TRUE)
    asreml::asreml.options(ai.sing = TRUE)
  }

  fit <- asreml::asreml(
    fixed = fixed,
    residual = resid,
    data = wide,
    na.action = asreml::na.method(y = "include", x = "include"),
    maxit = maxit,
    trace = FALSE,
    ...
  )
  if (!isTRUE(fit$converge)) {
    warning("asreml did not converge; treat the contrasts with caution.",
            call. = FALSE)
  }
  fit
}
