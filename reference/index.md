# Package index

## Designing the trial

Lay out the plots: strip or stacked, randomised or systematic.

- [`make_trial_design()`](https://www.zcao.space/ofeIntegrateR/reference/make_trial_design.md)
  : Lay out an on-farm strip trial

## Planning the sample

How many point samples the target precision needs, and where they go.

- [`kriging_sample_interval()`](https://www.zcao.space/ofeIntegrateR/reference/kriging_sample_interval.md)
  : How far apart should the samples be, and how many are needed?
- [`place_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/place_point_samples.md)
  : Choose where to put point samples

## Preparing the data

Getting an irregular yield-monitor cloud onto the lattice a spatial
model needs, and interpolating the sparse layer onto the same grid.

- [`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md)
  : Aggregate an irregular dense layer onto a regular lattice
- [`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md)
  : Krige sparse point-source samples onto a target grid

## Checking the point layer

Whether the samples are dense enough to reconstruct their own surface.

- [`cv_krige_surface()`](https://www.zcao.space/ofeIntegrateR/reference/cv_krige_surface.md)
  : Cross-validate a kriged point-source surface

## Fitting the spatial model

A REML engine with ASReml-style residual structures, written in base R,
plus the tests and predicted means that go with it.

- [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
  : Fit a spatial mixed model without asreml

- [`ofe_control()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_control.md)
  :

  Control parameters for
  [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)

- [`summary(`*`<ofe_fit>`*`)`](https://www.zcao.space/ofeIntegrateR/reference/summary.ofe_fit.md)
  :

  Summarise a fitted `ofe_fit`

- [`wald_tests()`](https://www.zcao.space/ofeIntegrateR/reference/wald_tests.md)
  : Wald tests of the fixed-effect terms

- [`ofe_means()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_means.md)
  : Predicted means for a fixed term

## Reporting the means

Predicted means with least significant differences and the a/b/c letters
a trial report prints beside them.

- [`ofe_lsd()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_lsd.md)
  : Predicted means with LSD letters
- [`compact_letters()`](https://www.zcao.space/ofeIntegrateR/reference/compact_letters.md)
  : Compact letter display from any set of pairwise comparisons

## Finding pseudo-environments

Cutting the trial into zones the spatial covariance can actually
support, and writing the residual structure that follows from them.

- [`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
  : Derive pseudo-environments from the spatial pattern of the dense
  layer
- [`adaptive_residual()`](https://www.zcao.space/ofeIntegrateR/reference/adaptive_residual.md)
  : Build a residual structure that matches the realised geometry

## Estimating the treatment effect

The two integration strategies, the free baseline they are judged
against, and a model-agnostic way to read the contrasts out.

- [`compare_integration()`](https://www.zcao.space/ofeIntegrateR/reference/compare_integration.md)
  : Compare an integrated analysis against the free baseline
- [`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)
  : Fit a treatment model using a kriged point-source covariate
- [`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md)
  : Fit a joint bivariate spatial model for dense and point-source data
- [`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
  : Extract a tidy table of fixed-effect estimates

## Simulating trials

A tidy lattice for quick tests, and the awkward shape a real harvester
produces for testing the full pipeline.

- [`simulate_ofe_trial()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_ofe_trial.md)
  : Simulate a synthetic OFE strip trial with a sparse point-source
  covariate
- [`simulate_yield_monitor()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_yield_monitor.md)
  : Simulate a trial as a yield monitor actually records it

## Package

- [`ofeIntegrateR`](https://www.zcao.space/ofeIntegrateR/reference/ofeIntegrateR-package.md)
  [`ofeIntegrateR-package`](https://www.zcao.space/ofeIntegrateR/reference/ofeIntegrateR-package.md)
  : ofeIntegrateR: Design, Sample and Analyse On-Farm Strip Experiments
