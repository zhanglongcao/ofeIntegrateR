# ofeIntegrateR 0.1.0

First release.

## New features

* `place_point_samples()` chooses where to put point samples under a named
  design — random, systematic grid, stratified or nested. The package Description
  promised sampling-design tools; this delivers them. Nested placement is a
  variogram-reconnaissance design and is deliberately clustered; grid and
  stratified placement are the ones to use when the samples will build a
  covariate surface.

* `cv_krige_surface()` cross-validates the kriged point-source surface, so a
  trial team can tell whether its point layer is dense enough to reconstruct
  itself. It is a sampling diagnostic, not a verdict on an analysis.

* `compare_integration()` fits the baseline and the integrated model on the same
  data and reports the treatment contrasts side by side, making the comparison
  the method should be judged on a single call.

* `simulate_yield_monitor()` generates a trial the way a harvester records it —
  GPS-referenced points along passes, position error, a clipped paddock corner,
  occasional missing passes, and point samples at arbitrary locations — so a
  workflow can be tested end to end through `grid_dense_layer()`.

* `grid_dense_layer()` aggregates an irregular dense layer — yield-monitor or
  proximal-sensor points along machinery passes — onto the complete rectangular
  lattice that a separable AR1 residual requires. Empty cells are retained as
  missing values rather than dropped, and `treat_purity` flags cells that
  straddle a treatment boundary.

* `fit_integrated_kriged(covariate = NULL)` fits the dense-layer-only baseline
  directly. This is the estimate an integrated analysis must beat, so it is now
  a first-class call rather than something to assemble by hand.

* `fit_integrated_kriged()` accepts several kriged covariates at once, so more
  than one point variable (soil nitrogen and phosphorus, say) can be integrated
  in a single model.

* New vignette, `vignette("ofe-integration")`: an end-to-end walkthrough from
  raw yield-monitor points to a treatment contrast, using open-source engines
  throughout.

## Infrastructure

* `sommer` moved from Imports to Suggests. It is only needed for the
  open-source joint-model engine and compiles C++ at install time, so users who
  only want the kriged-covariate route no longer pay for a build they will never
  call. `engine = "sommer"` now fails with a clear message if it is absent.

* Added a pkgdown site (`_pkgdown.yml`) with the reference index organised by
  pipeline stage, and a workflow to publish it to GitHub Pages.

* The R-CMD-check workflow now installs `knitr`, `rmarkdown` and `sommer`
  explicitly, so the vignette builds and the joint-model tests run in CI.

## Bug fixes

* `fit_integrated_kriged()` silently dropped all but the first covariate when
  passed more than one, because the formula was built with a vectorised
  `paste()` and no `collapse`. The extra covariates are now included.

* The `gls` engine stopped with "missing values in object" on gridded trial
  data, which routinely contains empty and excluded cells. Incomplete rows are
  now omitted unless the caller supplies its own `na.action`. (The `asreml`
  engine still keeps `NA` rows, as its separable AR1 residual needs the
  complete lattice.)

* Fixed the package URL and `BugReports` fields, which pointed at a GitHub
  account that no longer exists.
