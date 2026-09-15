# Changelog

## ofeIntegrateR 0.1.0

First release.

### New features

- [`place_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/place_point_samples.md)
  chooses where to put point samples under a named design — random,
  systematic grid, stratified or nested. The package Description
  promised sampling-design tools; this delivers them. Nested placement
  is a variogram-reconnaissance design and is deliberately clustered;
  grid and stratified placement are the ones to use when the samples
  will build a covariate surface.

- [`kriging_sample_interval()`](https://www.zcao.space/ofeIntegrateR/reference/kriging_sample_interval.md)
  answers the question that has to be settled before a sampling budget
  is set: given a variogram and a target precision, how far apart should
  the cores be, and how many does that imply for the trial’s area? It
  also reports the floor the nugget imposes — the precision no sampling
  density can beat — which is often the more useful number, because it
  says whether the target is worth budgeting for at all.

- [`make_trial_grid()`](https://www.zcao.space/ofeIntegrateR/reference/make_trial_grid.md)
  lays out a planned trial as a lattice with treatment strips, so
  sampling locations can be chosen before any yield data exist. Its
  output has the same shape as
  [`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md),
  so the same code works before and after harvest.

- [`cv_krige_surface()`](https://www.zcao.space/ofeIntegrateR/reference/cv_krige_surface.md)
  cross-validates the kriged point-source surface, so a trial team can
  tell whether its point layer is dense enough to reconstruct itself. It
  is a sampling diagnostic, not a verdict on an analysis.

- [`compare_integration()`](https://www.zcao.space/ofeIntegrateR/reference/compare_integration.md)
  fits the baseline and the integrated model on the same data and
  reports the treatment contrasts side by side, making the comparison
  the method should be judged on a single call.

- [`simulate_yield_monitor()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_yield_monitor.md)
  generates a trial the way a harvester records it — GPS-referenced
  points along passes, position error, a clipped paddock corner,
  occasional missing passes, and point samples at arbitrary locations —
  so a workflow can be tested end to end through
  [`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md).

- [`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md)
  aggregates an irregular dense layer — yield-monitor or proximal-sensor
  points along machinery passes — onto the complete rectangular lattice
  that a separable AR1 residual requires. Empty cells are retained as
  missing values rather than dropped, and `treat_purity` flags cells
  that straddle a treatment boundary.

- `fit_integrated_kriged(covariate = NULL)` fits the dense-layer-only
  baseline directly. This is the estimate an integrated analysis must
  beat, so it is now a first-class call rather than something to
  assemble by hand.

- [`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)
  accepts several kriged covariates at once, so more than one point
  variable (soil nitrogen and phosphorus, say) can be integrated in a
  single model.

- New vignette,
  [`vignette("ofe-integration")`](https://www.zcao.space/ofeIntegrateR/articles/ofe-integration.md):
  an end-to-end walkthrough from raw yield-monitor points to a treatment
  contrast, using open-source engines throughout.

### Infrastructure

- `sommer` moved from Imports to Suggests. It is only needed for the
  open-source joint-model engine and compiles C++ at install time, so
  users who only want the kriged-covariate route no longer pay for a
  build they will never call. `engine = "sommer"` now fails with a clear
  message if it is absent.

- Added a pkgdown site (`_pkgdown.yml`) with the reference index
  organised by pipeline stage, and a workflow to publish it to GitHub
  Pages.

- New article,
  [`vignette("real-data")`](https://www.zcao.space/ofeIntegrateR/articles/real-data.md):
  the full workflow on
  [`agridat::lasrosas.corn`](https://kwstat.github.io/agridat/reference/lasrosas.corn.html),
  an on-farm nitrogen experiment recorded by a yield monitor. It shows
  coordinate projection, choosing a cell size from `n_obs` and
  `treat_purity` rather than by convention, and a case where integration
  correctly turns out not to help — the point of reporting the baseline
  alongside it. `agridat` is a Suggests, and the article degrades
  gracefully when it is absent.

- The R-CMD-check workflow now installs `knitr`, `rmarkdown` and
  `sommer` explicitly, so the vignette builds and the joint-model tests
  run in CI.

- New `engine = "lme"` for
  [`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md),
  and `random` is now a first-class argument rather than something
  passed through `...`. This closes a gap that undermined the
  open-source path: `asreml` could fit a replicate or block effect,
  `gls` could not, and `lm` accepted the argument and silently discarded
  it. `lme` fits random effects together with the spatial correlation
  structure, so the open-source route is a genuine counterpart to the
  asreml one. `gls` and `lm` now raise an informative error instead of
  failing obscurely or ignoring the request. The asreml spelling
  `random = ~ rep` is translated to the grouping formula `lme` expects,
  so the same call works on both engines.

### Bug fixes

- [`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)
  silently dropped all but the first covariate when passed more than
  one, because the formula was built with a vectorised
  [`paste()`](https://rdrr.io/r/base/paste.html) and no `collapse`. The
  extra covariates are now included.

- The `gls` engine stopped with “missing values in object” on gridded
  trial data, which routinely contains empty and excluded cells.
  Incomplete rows are now omitted unless the caller supplies its own
  `na.action`. (The `asreml` engine still keeps `NA` rows, as its
  separable AR1 residual needs the complete lattice.)

- Fixed the package URL and `BugReports` fields, which pointed at a
  GitHub account that no longer exists.
