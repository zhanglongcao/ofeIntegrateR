#' Compare an integrated analysis against the free baseline
#'
#' Fits the treatment model twice on the same data --- once using the dense
#' layer alone, once with the kriged point-source covariate added --- and
#' reports the treatment contrasts side by side.
#'
#' This is the comparison the method should always be judged on. The dense layer
#' is already collected, so a spatial model of it costs nothing; point sampling
#' does. Reporting only the integrated fit hides how much of its accuracy came
#' from the covariate and how much was available for free. Report both.
#'
#' Note that a large, highly significant covariate coefficient is **not**
#' evidence that integration helped. With hundreds of grid cells the covariate
#' clears conventional significance almost automatically, including in
#' situations where it does nothing for the treatment estimate. Judge the
#' method by how the treatment contrasts and their standard errors move.
#'
#' @param data Data frame containing the gridded dense response, the treatment
#'   factor and the kriged covariate(s) --- typically the output of
#'   [krige_point_samples()].
#' @param response,treat Character; names of the response and treatment columns.
#' @param covariate Character vector; the kriged covariate column(s) to add in
#'   the integrated fit.
#' @param row,col Character; position columns used by the spatial residual.
#' @param engine Fitting engine, passed to [fit_integrated_kriged()].
#' @param ... Further arguments passed to [fit_integrated_kriged()].
#'
#' @return A data frame with one row per fixed-effect term, giving the estimate
#'   and standard error under each model and the change between them. The two
#'   fitted models are attached as `attr(, "models")`.
#'
#' @examples
#' sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 20, seed = 8)
#' kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
#' compare_integration(kr, response = "dense_response", treat = "treat",
#'                     covariate = "point_obs_kriged", engine = "lm")
#'
#' @export
compare_integration <- function(data,
                                 response,
                                 treat,
                                 covariate,
                                 row = "row",
                                 col = "col",
                                 engine = c("asreml", "gls", "lm"),
                                 ...) {
  engine <- match.arg(engine)
  if (length(covariate) == 0L) {
    stop("`covariate` must name at least one column; with none there is ",
         "nothing to compare against the baseline.", call. = FALSE)
  }

  baseline <- fit_integrated_kriged(data, response = response, treat = treat,
                                     covariate = NULL, row = row, col = col,
                                     engine = engine, ...)
  integrated <- fit_integrated_kriged(data, response = response, treat = treat,
                                       covariate = covariate, row = row,
                                       col = col, engine = engine, ...)

  b <- extract_fixed_effects(baseline)
  i <- extract_fixed_effects(integrated)
  names(b)[names(b) %in% c("estimate", "se")] <-
    paste0("baseline_", c("estimate", "se"))
  names(i)[names(i) %in% c("estimate", "se")] <-
    paste0("integrated_", c("estimate", "se"))

  out <- merge(b, i, by = "term", all = TRUE, sort = FALSE)
  out$estimate_change <- out$integrated_estimate - out$baseline_estimate
  # Negative = the integrated fit is more precise for this term.
  out$se_change <- out$integrated_se - out$baseline_se
  rownames(out) <- NULL

  attr(out, "models") <- list(baseline = baseline, integrated = integrated)
  out
}
