# ofeIntegrateR 0.1.0

First release.

## New features

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
