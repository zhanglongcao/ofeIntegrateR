#' Cross-validate a kriged point-source surface
#'
#' Leave-one-out cross-validation of the point-source layer: each sample is
#' dropped in turn, predicted from the rest, and compared with its observed
#' value. It answers whether the samples are dense enough, and well enough
#' behaved, to reconstruct their own surface.
#'
#' Treat this as a **sampling** diagnostic, not a verdict on an analysis. In the
#' simulation work behind this package a better-validated surface delivered a
#' larger average benefit from integration, but the probability of benefiting on
#' any individual trial barely moved. A low or negative `r2` is a reason to
#' collect more samples next season; it does not by itself condemn the current
#' model. Note also that cross-validation measures whether the surface is
#' *well estimated*, not whether it is *relevant to yield* — a soil property can
#' be mapped perfectly and still explain nothing about the treatment response.
#'
#' @param point_data Data frame of point-source observations.
#' @param value Character; name of the observed value column.
#' @param coords Character vector of length 2 naming the coordinate columns.
#' @param vgm_model Optional variogram model from [gstat::vgm()]. By default the
#'   same exponential model [krige_point_samples()] would fit is used, so the
#'   diagnostic matches the surface actually entering the analysis.
#' @param nmax Maximum neighbours per prediction, passed to [gstat::krige.cv()].
#'
#' @return A list with `rmse` (root mean squared prediction error), `r2`
#'   (proportion of variance explained, which is negative when kriging predicts
#'   worse than the sample mean), `n` (number of points) and `residuals`.
#'   The variogram used is attached as `attr(, "variogram")`.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 25, seed = 4)
#' cv <- cv_krige_surface(sim$point_samples, value = "point_obs")
#' cv$rmse
#' cv$r2
#'
#' @export
cv_krige_surface <- function(point_data,
                              value,
                              coords = c("col", "row"),
                              vgm_model = NULL,
                              nmax = 30) {
  if (!value %in% names(point_data)) {
    stop("`value` column '", value, "' not found in `point_data`.", call. = FALSE)
  }
  if (!all(coords %in% names(point_data))) {
    stop("`coords` columns must be present in `point_data`: ",
         paste(setdiff(coords, names(point_data)), collapse = ", "), call. = FALSE)
  }
  obs <- point_data[[value]]
  keep <- is.finite(obs)
  point_data <- point_data[keep, , drop = FALSE]
  obs <- obs[keep]
  n <- nrow(point_data)
  if (n < 4L) {
    stop("At least 4 point observations are needed to cross-validate.",
         call. = FALSE)
  }

  coord_formula <- stats::as.formula(paste("~", coords[1], "+", coords[2]))
  sp_pts <- point_data
  sp::coordinates(sp_pts) <- coord_formula
  value_formula <- stats::as.formula(paste(value, "~ 1"))

  if (is.null(vgm_model)) {
    extent <- max(diff(range(point_data[[coords[1]]])),
                  diff(range(point_data[[coords[2]]])))
    start <- gstat::vgm(psill = stats::var(obs), model = "Exp",
                        range = extent / 4, nugget = 0.1)
    emp <- gstat::variogram(value_formula, sp_pts)
    vgm_model <- tryCatch(gstat::fit.variogram(emp, start),
                          error = function(e) start,
                          warning = function(w) start)
  }

  cv <- gstat::krige.cv(value_formula, sp_pts, model = vgm_model,
                        nfold = n, nmax = nmax, verbose = FALSE)

  res <- cv$residual
  sst <- sum((obs - mean(obs))^2)
  out <- list(
    rmse = sqrt(mean(res^2, na.rm = TRUE)),
    r2   = if (sst > 0) 1 - sum(res^2, na.rm = TRUE) / sst else NA_real_,
    n    = n,
    residuals = res
  )
  attr(out, "variogram") <- vgm_model
  out
}
