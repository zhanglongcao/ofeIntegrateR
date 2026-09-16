#' ofeIntegrateR: Design, Sample and Analyse On-Farm Strip Experiments
#'
#' A workflow for on-farm experimentation (OFE) strip trials that runs
#' end to end without a commercial licence.
#'
#' @section Design the trial:
#' [make_trial_design()] lays out treatment plots in replicate blocks, as a
#' strip layout or two stacked tiers with a buffer, with randomised or
#' systematic treatment order.
#'
#' @section Plan the sampling:
#' [kriging_sample_interval()] turns a variogram and a target precision into a
#' sampling interval and a core count -- and reports the precision floor the
#' nugget imposes. [place_point_samples()] chooses the locations.
#'
#' @section Prepare the data:
#' [grid_dense_layer()] snaps an irregular yield-monitor cloud onto the
#' complete lattice a separable residual needs. [krige_point_samples()]
#' interpolates the sparse layer onto the same grid, and [cv_krige_surface()]
#' says whether it was dense enough to be worth it.
#'
#' @section Analyse:
#' [fit_ofe()] fits a spatial mixed model by REML with ASReml-style residual
#' structures -- `ar1()`, `id()`, `exp()`, `diag()`, and `dsum()` sections --
#' in base R. [wald_tests()] and [ofe_means()] are the analogues of
#' `wald.asreml()` and `predict.asreml()`. On the structures both support it
#' agrees with `asreml::asreml()` to several significant figures.
#'
#' @section Report the means:
#' [ofe_lsd()] gives the predicted means with their least significant
#' difference and the a/b/c letters, matching `agricolae::LSD.test()` and
#' `HSD.test()` on a balanced design. [compact_letters()] letters pairwise
#' comparisons from any other source.
#'
#' @section Pseudo-environments:
#' Two ways to find them, for two different questions.
#' [partition_pseudo_env()] slices a strip trial across its length from the
#' pattern in the dense layer, using the fitted correlation range to stop the
#' partition running past the scale the data can support; every treatment stays
#' in every zone. [partition_paddock()] clusters environmental covariates --
#' elevation, EM38, soil tests -- into contiguous regions of any shape, using a
#' minimum spanning tree so that a zone cannot come back as scattered cells the
#' way k-means leaves it. [adaptive_residual()] writes the matching `dsum()`
#' residual formula for either, demoting `ar1()` to `id()` in any zone that
#' lacks the extent to support it.
#'
#' @section Integrate a point-source layer:
#' Two strategies, each with an open-source path:
#'
#' \describe{
#'   \item{Kriged-covariate ([fit_integrated_kriged()])}{Krige the sparse
#'     variable onto the trial grid, then include the surface as a fixed
#'     covariate in a spatial mixed model of the dense response. Pass
#'     `covariate = NULL` for the dense-layer-only baseline the integrated fit
#'     has to beat.}
#'   \item{Joint bivariate model ([fit_integrated_joint()])}{Fit the dense
#'     response and the point-source measurements jointly as two stacked
#'     layers sharing a spatial random effect, with no separate kriging step.}
#' }
#'
#' [compare_integration()] fits the baseline and the integrated model together,
#' and [extract_fixed_effects()] reads contrasts out of any of the fits.
#'
#' @section Simulate:
#' [simulate_ofe_trial()] for a tidy lattice, [simulate_yield_monitor()] for
#' the awkward shape a real harvester produces.
#'
#' Developed for the AAGI-CU-RD-OFE GRDC project
#' (\dQuote{Development of processes to integrate point-source data and
#' high-resolution data}).
#'
#' @keywords internal
#' @importFrom stats coef fitted logLik nobs residuals vcov
"_PACKAGE"
