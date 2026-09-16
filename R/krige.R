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
#'   starting values for [gstat::fit.variogram()]. Defaults to the constrained
#'   fit from [ofe_variogram()], which cannot fail to converge: given the range
#'   it solves for nugget and sill in closed form and searches the range on a
#'   grid. `gstat::fit.variogram()` then refines it where it can; where it
#'   cannot, that fit is used unchanged and `attr(, "variogram_source")` on the
#'   result says so.
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
    # Start from the package's own constrained fit. The previous defaults were
    # a fixed nugget of 0.1 whatever the variable's scale, which for anything
    # measured in hundreds (brightness, EC) is effectively zero, and gstat's
    # weighted least squares would wander for 200 iterations and give up.
    # ofe_variogram() cannot fail to converge: given the range it solves for
    # the nugget and sill in closed form and searches the range on a grid.
    v0 <- tryCatch(
      ofe_variogram(point_data, value = value, x = coords[1], y = coords[2]),
      error = function(e) NULL)
    vgm_start <- if (!is.null(v0)) {
      gstat::vgm(psill = v0$psill, nugget = v0$nugget, range = v0$range,
                 model = switch(v0$model, exponential = "Exp",
                                spherical = "Sph", gaussian = "Gau", "Exp"))
    } else {
      # Too few samples for a variogram at all; fall back to scale-aware
      # heuristics rather than to a constant.
      vv <- stats::var(point_data[[value]], na.rm = TRUE)
      extent <- max(diff(range(point_data[[coords[1]]])),
                    diff(range(point_data[[coords[2]]])))
      gstat::vgm(psill = 0.9 * vv, model = "Exp", range = extent / 4,
                 nugget = 0.1 * vv)
    }
  }

  coord_formula <- stats::as.formula(paste("~", coords[1], "+", coords[2]))
  point_sp <- point_data
  sp::coordinates(point_sp) <- coord_formula

  value_formula <- stats::as.formula(paste(value, "~ 1"))
  vgm_emp <- gstat::variogram(value_formula, point_sp)
  # gstat's fit refines the starting values; when it cannot, the starting
  # values are themselves a fitted model now, so falling back to them is a
  # substitution rather than a degradation and does not warrant a warning.
  # Muffling the warning is not enough: fit.variogram still returns its
  # non-converged iterate, so the substitution has to be made explicitly or
  # the message would claim a fallback that never happened.
  converged <- TRUE
  vgm_fit <- withCallingHandlers(
    tryCatch(
      gstat::fit.variogram(vgm_emp, vgm_start),
      error = function(e) {
        converged <<- FALSE
        vgm_start
      }
    ),
    warning = function(w) {
      if (grepl("No convergence", conditionMessage(w), fixed = TRUE)) {
        converged <<- FALSE
        invokeRestart("muffleWarning")
      }
    }
  )
  # Silent by design. The starting values are a fitted model, not a guess, so
  # gstat declining to refine them is an implementation detail rather than
  # something wrong -- and on small point sets it declines often enough that a
  # message on every call would be noise. Which fit was used is recorded on the
  # result instead, so it stays inspectable.
  refined <- converged && all(is.finite(vgm_fit$psill)) &&
    all(vgm_fit$psill >= 0)
  if (!refined) vgm_fit <- vgm_start

  newdata_sp <- newdata
  sp::coordinates(newdata_sp) <- coord_formula

  krige_out <- gstat::krige(value_formula, point_sp, newdata_sp,
                             model = vgm_fit, nmax = nmax, debug.level = 0)

  out <- newdata
  out[[paste0(value, "_kriged")]] <- krige_out$var1.pred
  out[[paste0(value, "_kriged_var")]] <- krige_out$var1.var
  attr(out, "variogram") <- vgm_fit
  attr(out, "variogram_source") <- if (refined) "gstat::fit.variogram" else
    "ofeIntegrateR constrained fit"
  out
}
