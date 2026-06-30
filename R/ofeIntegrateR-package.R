#' ofeIntegrateR: Integrate Point-Source and High-Resolution Spatial Data
#' in On-Farm Experiments
#'
#' Tools for integrating sparse point-source measurements (e.g. soil cores,
#' tissue samples, disease ratings) with dense spatial covariates (e.g.
#' EM38 surveys, yield maps) in on-farm experimentation (OFE) strip trials.
#'
#' Two integration strategies are provided:
#'
#' \describe{
#'   \item{Kriged-covariate ([fit_integrated_kriged()])}{Krige the sparse
#'     point-source variable onto the trial grid with [krige_point_samples()],
#'     then include the kriged surface as a fixed covariate in a spatial
#'     mixed model of the dense response.}
#'   \item{Joint bivariate model ([fit_integrated_joint()])}{Fit the dense
#'     response and the point-source measurements jointly as two stacked
#'     response layers sharing a spatial random effect, avoiding a separate
#'     kriging step.}
#' }
#'
#' [simulate_ofe_trial()] generates synthetic test data for both strategies,
#' and [extract_fixed_effects()] provides a model-agnostic way to retrieve
#' fixed-effect (e.g. treatment contrast) estimates from either an `lm` or
#' an `asreml` fit.
#'
#' Developed for Milestone 4 of the AAGI-CU-RD-OFE GRDC project
#' (\dQuote{Development of processes to integrate point-source data and
#' high-resolution data}).
#'
#' @keywords internal
"_PACKAGE"
