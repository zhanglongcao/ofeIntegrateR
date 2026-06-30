#' Krige sparse point-source samples onto a target grid
#'
#' Fits an exponential variogram to sparse point-source observations (e.g.
#' soil cores, tissue samples) and uses ordinary kriging to predict values
#' across a target grid (e.g. the full trial extent or a dense covariate
#' grid). This is the first step of the "kriged-covariate" integration
#' strategy implemented in [fit_integrated_kriged()].
#'
#' @param point_data Data frame of sparse point-source observations,
#'   containing the coordinate columns and the value column.
#' @param newdata Data frame giving the target locations to predict onto
#'   (e.g. every cell of the trial grid), containing the same coordinate
#'   columns as `point_data`.
#' @param value Character; name of the column in `point_data` holding the
#'   observed point-source value.
#' @param coords Character vector of length 2 giving the coordinate column
#'   names present in both `point_data` and `newdata`.
#' @param vgm_start A variogram model object from [gstat::vgm()] giving
#'   starting values for [gstat::fit.variogram()]. Defaults to an
#'   exponential model with psill = 1, range = `max(spatial extent) / 4`,
#'   nugget = 0.1.
#' @param nmax Maximum number of nearest observations used for each
#'   kriging prediction (passed to [gstat::krige()]).
#'
#' @return `newdata` with two columns appended: `<value>_kriged` (the
#'   kriging prediction) and `<value>_kriged_var` (the kriging variance).
#'   The fitted variogram model is attached as `attr(., "variogram")`.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 12, seed = 1)
#' krieged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
#' head(krieged)
#' attr(krieged, "variogram")
#'
#' @export
krige_point_samples <- function(point_data,
                                 newdata,
                                 value,
                                 coords = c("col", "row"),
                                 vgm_start = NULL,
                                 nmax = 30) {
  if (!value %in% names(point_data)) {
    stop("`value` column '", value, "' not found in `point_data`.", call. = FALSE)
  }
  if (!all(coords %in% names(point_data)) || !all(coords %in% names(newdata))) {
    stop("`coords` columns must be present in both `point_data` and `newdata`.",
         call. = FALSE)
  }

  if (is.null(vgm_start)) {
    extent <- max(diff(range(point_data[[coords[1]]])),
                  diff(range(point_data[[coords[2]]])))
    vgm_start <- gstat::vgm(psill = stats::var(point_data[[value]], na.rm = TRUE),
                             model = "Exp", range = extent / 4, nugget = 0.1)
  }

  coord_formula <- stats::as.formula(paste("~", coords[1], "+", coords[2]))
  point_sp <- point_data
  sp::coordinates(point_sp) <- coord_formula

  value_formula <- stats::as.formula(paste(value, "~ 1"))
  vgm_emp <- gstat::variogram(value_formula, point_sp)
  vgm_fit <- tryCatch(
    gstat::fit.variogram(vgm_emp, vgm_start),
    error = function(e) {
      warning("Variogram fitting failed (", conditionMessage(e),
              "); falling back to starting values.", call. = FALSE)
      vgm_start
    }
  )

  newdata_sp <- newdata
  sp::coordinates(newdata_sp) <- coord_formula

  krige_out <- gstat::krige(value_formula, point_sp, newdata_sp,
                             model = vgm_fit, nmax = nmax, debug.level = 0)

  out <- newdata
  out[[paste0(value, "_kriged")]] <- krige_out$var1.pred
  out[[paste0(value, "_kriged_var")]] <- krige_out$var1.var
  attr(out, "variogram") <- vgm_fit
  out
}
