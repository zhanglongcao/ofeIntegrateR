# ofeIntegrateR: Design, Sample and Analyse On-Farm Strip Experiments

A workflow for on-farm experimentation (OFE) strip trials that runs end
to end without a commercial licence.

## Design the trial

[`make_trial_design()`](https://www.zcao.space/ofeIntegrateR/reference/make_trial_design.md)
lays out treatment plots in replicate blocks, as a strip layout or two
stacked tiers with a buffer, with randomised or systematic treatment
order.

## Plan the sampling

[`kriging_sample_interval()`](https://www.zcao.space/ofeIntegrateR/reference/kriging_sample_interval.md)
turns a variogram and a target precision into a sampling interval and a
core count – and reports the precision floor the nugget imposes.
[`place_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/place_point_samples.md)
chooses the locations.

## Look at it

[`ofe_map()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_map.md)
draws any column over the trial – yield, a kriged soil surface,
elevation, the zones, the residuals – and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods cover
the trial design
([`plot.ofe_design()`](https://www.zcao.space/ofeIntegrateR/reference/plot.ofe_design.md)),
the zones
([`plot.ofe_zones()`](https://www.zcao.space/ofeIntegrateR/reference/plot.ofe_zones.md)),
the variogram
([`plot.ofe_variogram()`](https://www.zcao.space/ofeIntegrateR/reference/plot.ofe_variogram.md))
and a fitted model's diagnostics
([`plot.ofe_fit()`](https://www.zcao.space/ofeIntegrateR/reference/plot.ofe_fit.md)).
[`ofe_palette()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_palette.md)
holds the colours, chosen against a colour-vision validator rather than
by eye.

## Prepare the data

[`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md)
snaps an irregular yield-monitor cloud onto the complete lattice a
separable residual needs.
[`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md)
interpolates the sparse layer onto the same grid, and
[`cv_krige_surface()`](https://www.zcao.space/ofeIntegrateR/reference/cv_krige_surface.md)
says whether it was dense enough to be worth it.
[`ofe_variogram()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_variogram.md)
fits and draws the variogram those steps depend on: the range says how
far one core speaks for, and the nugget says how much variation no
sampling density will ever resolve.

## Analyse

[`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
fits a spatial mixed model by REML with ASReml-style residual structures
– `ar1()`, `id()`, [`exp()`](https://rdrr.io/r/base/Log.html),
[`diag()`](https://rdrr.io/r/base/diag.html), and `dsum()` sections – in
base R.
[`wald_tests()`](https://www.zcao.space/ofeIntegrateR/reference/wald_tests.md)
and
[`ofe_means()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_means.md)
are the analogues of `wald.asreml()` and `predict.asreml()`. On the
structures both support it agrees with `asreml::asreml()` to several
significant figures.

## Report the means

[`ofe_lsd()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_lsd.md)
gives the predicted means with their least significant difference and
the a/b/c letters, matching `agricolae::LSD.test()` and `HSD.test()` on
a balanced design.
[`compact_letters()`](https://www.zcao.space/ofeIntegrateR/reference/compact_letters.md)
letters pairwise comparisons from any other source.

## Pseudo-environments

Two ways to find them, for two different questions.
[`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
slices a strip trial across its length from the pattern in the dense
layer, using the fitted correlation range to stop the partition running
past the scale the data can support; every treatment stays in every
zone.
[`partition_paddock()`](https://www.zcao.space/ofeIntegrateR/reference/partition_paddock.md)
zones a paddock on environmental covariates – elevation, EM38, soil
tests – into rectangular blocks by default, so the zones are something a
machine can drive and a sampling grid can follow, with contiguous
free-form zones available where the soil matters more than the shape.
[`adaptive_residual()`](https://www.zcao.space/ofeIntegrateR/reference/adaptive_residual.md)
writes the matching `dsum()` residual formula for either, demoting
`ar1()` to `id()` in any zone that lacks the extent to support it.

## Integrate a point-source layer

Two strategies, each with an open-source path:

- Kriged-covariate
  ([`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)):

  Krige the sparse variable onto the trial grid, then include the
  surface as a fixed covariate in a spatial mixed model of the dense
  response. Pass `covariate = NULL` for the dense-layer-only baseline
  the integrated fit has to beat.

- Joint bivariate model
  ([`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md)):

  Fit the dense response and the point-source measurements jointly as
  two stacked layers sharing a spatial random effect, with no separate
  kriging step.

[`compare_integration()`](https://www.zcao.space/ofeIntegrateR/reference/compare_integration.md)
fits the baseline and the integrated model together, and
[`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
reads contrasts out of any of the fits.

## Simulate

[`simulate_paddock()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_paddock.md)
builds a paddock whose truth is known: a correlated yield-potential
surface, covariate layers related to it by a chosen correlation,
pseudo-environments with their own treatment response, and a trial laid
into it – so the zoning and the analysis can be scored rather than
merely run.
[`simulate_ofe_trial()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_ofe_trial.md)
gives a tidy lattice for quick tests, and
[`simulate_yield_monitor()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_yield_monitor.md)
the awkward shape a real harvester produces.

Developed for the AAGI-CU-RD-OFE GRDC project (“Development of processes
to integrate point-source data and high-resolution data”).

## See also

Useful links:

- <https://github.com/zhanglongcao/ofeIntegrateR>

- <https://www.zcao.space/ofeIntegrateR/>

- Report bugs at <https://github.com/zhanglongcao/ofeIntegrateR/issues>

## Author

**Maintainer**: Zhanglong Cao <zhanglong.cao@curtin.edu.au>
