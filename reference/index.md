# Package index

## Planning the trial and the sample

Before the season: lay out the trial, work out how many cores the target
precision needs, and choose where they go.

- [`make_trial_grid()`](https://www.zcao.space/ofeIntegrateR/reference/make_trial_grid.md)
  : Build a candidate grid for a planned trial
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
  : ofeIntegrateR: Integrate Point-Source and High-Resolution Spatial
  Data in On-Farm Experiments
