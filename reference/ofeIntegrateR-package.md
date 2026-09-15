# ofeIntegrateR: Integrate Point-Source and High-Resolution Spatial Data in On-Farm Experiments

Tools for integrating sparse point-source measurements (e.g. soil cores,
tissue samples, disease ratings) with dense spatial covariates (e.g.
EM38 surveys, yield maps) in on-farm experimentation (OFE) strip trials.

## Details

Two integration strategies are provided:

- Kriged-covariate
  ([`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)):

  Krige the sparse point-source variable onto the trial grid with
  [`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md),
  then include the kriged surface as a fixed covariate in a spatial
  mixed model of the dense response.

- Joint bivariate model
  ([`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md)):

  Fit the dense response and the point-source measurements jointly as
  two stacked response layers sharing a spatial random effect, avoiding
  a separate kriging step.

[`simulate_ofe_trial()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_ofe_trial.md)
generates synthetic test data for both strategies, and
[`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
provides a model-agnostic way to retrieve fixed-effect (e.g. treatment
contrast) estimates from either an `lm` or an `asreml` fit.

Developed for Milestone 4 of the AAGI-CU-RD-OFE GRDC project
(“Development of processes to integrate point-source data and
high-resolution data”).

## See also

Useful links:

- <https://github.com/zhanglongcao/ofeIntegrateR>

- <https://www.zcao.space/ofeIntegrateR/>

- Report bugs at <https://github.com/zhanglongcao/ofeIntegrateR/issues>

## Author

**Maintainer**: Zhanglong Cao <zhanglong.cao@curtin.edu.au>
