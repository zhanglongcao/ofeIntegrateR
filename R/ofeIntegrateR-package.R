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
#' @section Look at it:
#' [ofe_map()] draws any column over the trial -- yield, a kriged soil surface,
#' elevation, the zones, the residuals -- and `plot()` methods cover the trial
#' design ([plot.ofe_design()]), the zones ([plot.ofe_zones()]), the variogram
#' ([plot.ofe_variogram()]) and a fitted model's diagnostics
#' ([plot.ofe_fit()]). [ofe_palette()] holds the colours, chosen against a
#' colour-vision validator rather than by eye.
#'
#' @section Prepare the data:
#' [grid_dense_layer()] snaps an irregular yield-monitor cloud onto the
#' complete lattice a separable residual needs. [krige_point_samples()]
#' interpolates the sparse layer onto the same grid, and [cv_krige_surface()]
#' says whether it was dense enough to be worth it. [ofe_variogram()] fits and
#' draws the variogram those steps depend on: the range says how far one core
#' speaks for, and the nugget says how much variation no sampling density will
#' ever resolve.
#'
#' @section Analyse:
#' [fit_ofe()] fits a spatial mixed model by REML with separable residual
#' structures -- `ar1()`, `id()`, `exp()`, `diag()`, and `dsum()` sections --
#' in base R. [wald_tests()] tests each fixed term as a whole and [ofe_means()]
#' gives predicted means and contrasts. The residual spellings follow ASReml-R,
#' and on the structures both support the estimates agree to several
#' significant figures, so a script written against either runs on the other.
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
#' in every zone. [partition_paddock()] zones a paddock on environmental
#' covariates -- elevation, EM38, soil tests -- into rectangular blocks by
#' default, so the zones are something a machine can drive and a sampling grid
#' can follow, with contiguous free-form zones available where the soil matters
#' more than the shape. [adaptive_residual()] writes the matching `dsum()`
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
#' [simulate_paddock()] builds a paddock whose truth is known: a correlated
#' yield-potential surface, covariate layers related to it by a chosen
#' correlation, pseudo-environments with their own treatment response, and a
#' trial laid into it -- so the zoning and the analysis can be scored rather
#' than merely run. [simulate_ofe_trial()] gives a tidy lattice for quick
#' tests, and [simulate_yield_monitor()] the awkward shape a real harvester
#' produces.
#'
#' Developed for the AAGI-CU-RD-OFE GRDC project
#' (\dQuote{Development of processes to integrate point-source data and
#' high-resolution data}).
#'
#' @keywords internal
#' @importFrom stats coef fitted logLik nobs residuals vcov
"_PACKAGE"
