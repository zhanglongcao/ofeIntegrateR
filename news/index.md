# Changelog

## ofeIntegrateR 0.1.0

First release.

### New features

- [`simulate_paddock()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_paddock.md)
  builds a paddock whose truth is known, which the package had no way to
  do: the existing simulators produced a trial but never covariate
  layers (so
  [`partition_paddock()`](https://www.zcao.space/ofeIntegrateR/reference/partition_paddock.md)
  had nothing to zone), never zones (so
  [`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
  had no answer to be checked against), and never separated yield
  potential from treatment response (so nothing could be attributed). It
  draws a correlated yield-potential surface, then each named covariate
  as a mixture of that surface and its own structure in a proportion you
  set – a covariate unrelated to yield is not worth zoning on, and one
  identical to it is an unrealistically easy test. Pseudo-environments
  cut across the strips and scale the treatment response zone by zone,
  so a zone-by-treatment interaction is genuinely present or genuinely
  absent, and both cases can be tested. The settings and the realised
  component variances are attached as `attr(, "truth")`.

- [`ofe_variogram()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_variogram.md)
  fits and draws the variogram, which the package had been kriging from
  without ever showing. It bins every pair of samples by distance, fits
  exponential, spherical and Gaussian models by a constrained search,
  and reports the nugget, sill and practical range with the usual
  classification of spatial dependence. Two things it does that a plain
  variogram call does not: `trend = ~ treat` strips a mean effect first,
  since a treatment step keeps the variogram climbing and inflates the
  range; and when the fitted practical range is longer than the largest
  lag the samples cover, it says so, in
  [`print()`](https://rdrr.io/r/base/print.html) and on the plot, rather
  than reporting a number the data never supported.
  [`kriging_sample_interval()`](https://www.zcao.space/ofeIntegrateR/reference/kriging_sample_interval.md)
  now takes the fitted object directly.

- [`ofe_map()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_map.md)
  draws any column over the trial, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods cover
  the design, the zones, the variogram and a fitted model.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on an
  `ofe_fit` gives the diagnostic the package most needed: a spatial
  model is fitted to absorb the field’s pattern, so the question is
  whether any is left, and the residual map and residual variogram
  answer it.

- [`ofe_palette()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_palette.md)
  holds the plotting colours. They were chosen by running candidate sets
  through a colour-vision validator against the surface an R device
  draws on, not by eye: the categorical order is pairwise-separable up
  to six treatments, the ordinal zone ramp up to seven zones, and label
  ink is picked per fill by whichever of dark or light actually has more
  contrast.
  [`ofe_map()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_map.md)
  will not select a diverging scale on its own – a variable that merely
  contains negative values is not a signed one, and a two-hue scale on
  it invents a midpoint the data does not have.

- [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
  fits a linear mixed model by residual maximum likelihood with a
  separable, ASReml-style residual structure, written in base R. It
  exists because the analysis was the one step of the pipeline that
  still needed a commercial licence: the model at the end of an OFE
  workflow is `yield ~ treat` with an `ar1(row):ar1(col)` residual, and
  until now that meant asreml or an approximation. `residual` takes the
  asreml spellings – `id()`, `ar1()`,
  [`exp()`](https://rdrr.io/r/base/Log.html),
  [`diag()`](https://rdrr.io/r/base/diag.html), and
  `dsum(~ struct | s, levels = )` sections – combined with `:`, so a
  script written against asreml runs unchanged. On the structures both
  support the two agree to several significant figures, on variance
  components, fixed effects and their standard errors alike.
  [`wald_tests()`](https://www.zcao.space/ofeIntegrateR/reference/wald_tests.md)
  and
  [`ofe_means()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_means.md)
  are the analogues of `wald.asreml()` and `predict.asreml()`, the
  latter also returning pairwise contrasts with their standard errors of
  difference.

  Unlike asreml, cells with a missing response are dropped rather than
  padded: the correlation is evaluated from the positions of the
  observations that remain, and a row absent from the data still
  contributes its lag instead of being closed up. The engine is dense –
  it factorises an n-by-n matrix at every iteration – so it is meant for
  OFE-sized problems, and refuses more than `ofe_control(max_n = )`
  observations rather than appearing to hang.

- [`ofe_lsd()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_lsd.md)
  reports the predicted means with their least significant difference
  and a compact letter display – the a/b/c annotation a trial report
  prints beside its treatment means. On a balanced design it reproduces
  `agricolae::LSD.test()` exactly, LSD value and letters alike, and
  `adjust = "tukey"` reproduces `agricolae::HSD.test()`. Its default
  differs from both on the data this package is for: each pair is judged
  on its own standard error of difference, because under a spatial model
  neighbouring strips are compared more precisely than distant ones, and
  a single average LSD hides that. `use = "lsd"` gives the classic
  single LSD for a published table, and the value is returned alongside
  so it can be quoted. `by = "zone"` letters the means within each
  pseudo-environment, with comparisons and LSDs kept inside a zone.

- [`compact_letters()`](https://www.zcao.space/ofeIntegrateR/reference/compact_letters.md)
  turns any set of pairwise comparisons into letters, so output from
  `asreml::predict()`, `emmeans` or a hand-typed table can be annotated
  the same way. It takes a p-value matrix, a `contrast` column of the
  form `"B - A"`, or explicit `level1`/`level2` columns. Letters come
  from the maximal cliques of the not-different graph, which is the
  minimal correct display; a pair with no p-value raises a warning
  rather than passing silently as “different”.

- [`ofe_means()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_means.md)
  gains `by`, for means within each level of a second factor, and
  `adjust`, for multiplicity-adjusted pairwise p-values (`"tukey"`,
  `"sidak"`, or any \[stats::p.adjust()\] method). Its pairwise table
  now names the two levels of each comparison in `level1` and `level2`
  as well as in `contrast`.

- [`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
  derives pseudo-environments from the dense layer instead of taking
  them as given. It removes the treatment signal, collapses the residual
  field to a profile along the trial, and segments that profile
  optimally by dynamic programming. The spatial covariance is what stops
  it inventing zones: the practical range of a fitted exponential
  variogram becomes the minimum zone width and caps the number of zones
  at `floor(trial length / range)`, because a region narrower than the
  correlation range is one realisation of the same surface rather than a
  distinct environment. Within that cap, BIC chooses. The derivation is
  returned in full – the fitted range, the breaks, the BIC table and the
  profile – so the answer can be argued with rather than merely
  accepted.

- [`partition_paddock()`](https://www.zcao.space/ofeIntegrateR/reference/partition_paddock.md)
  zones a paddock into pseudo-environments from environmental covariates
  – elevation, EM38 or gamma survey, soil test grids – rather than from
  the yield being analysed. `covariates` names the columns of `data` to
  zone on and is required: there is no default, because which layers
  define a zone is agronomy rather than something the function can
  guess.

  Three methods, differing only in the shape a zone is allowed to be,
  which is a practical decision rather than a statistical one.
  `"rectangle"` (the default) cuts the paddock by lines running the full
  width or length of the region being split, each placed where it
  removes most within-zone variance; cutting a rectangle across leaves
  two rectangles, so every zone is a rectangle however many cuts are
  made, and the zone summary reports each block’s corners so it can be
  marked out. It is the two-dimensional version of what
  [`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
  does along one axis. `"skater"` gives contiguous zones of any shape,
  from a minimum spanning tree over the neighbourhood graph pruned one
  edge at a time (Assuncao et al. 2006) – right when a boundary really
  runs at an angle, at the cost of zones that send fingers out between
  their neighbours. `"kmeans"` is not contiguous at all and is kept only
  for comparison: it assigns each cell to the zone its soil resembles
  wherever the cell sits, so zones arrive as confetti. The `patches`
  column counts the connected pieces of each zone, so the difference can
  be seen rather than taken on trust.

  Being freer does not make skater fit better: both it and the rectangle
  method are greedy, and on the test paddock the rectangles explain more
  of the covariate variance at the same k (r2 0.78 against 0.73). The
  `r2` column of the returned table is there to be compared rather than
  assumed.

  With `k` unset, zones are added while each earns its keep: the default
  stops at the first `k` whose successor would explain less than
  `min_gain` (5%) more of the covariate variance, and a paddock uniform
  in its covariates comes back as one zone. The Calinski-Harabasz
  criterion is available, documented with the caveat that on a smoothly
  varying covariate it usually has no interior maximum and simply picks
  `k_max`. Passing `treat` reports each zone’s treatment coverage and
  warns where a zone cannot support a `zone:treat` term.

- [`adaptive_residual()`](https://www.zcao.space/ofeIntegrateR/reference/adaptive_residual.md)
  writes the `dsum()` residual formula a set of zones implies, giving
  each zone `ar1()` only in the dimensions where it has extent and
  `id()` elsewhere. Without this, a zone one row deep asks for an
  unidentifiable AR1 parameter and the fit fails. The formula is plain
  text that both
  [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
  and `asreml::asreml()` accept.

- [`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)
  and
  [`compare_integration()`](https://www.zcao.space/ofeIntegrateR/reference/compare_integration.md)
  now default to `engine = "ofe"`, so the documented workflow fits the
  AR1xAR1 model with no licence. `engine = "asreml"` still fits the
  identical model where licensed, and
  [`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
  reads contrasts out of either.

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

- [`make_trial_design()`](https://www.zcao.space/ofeIntegrateR/reference/make_trial_design.md)
  lays out the trial itself: treatment plots of a given width and length
  in replicate blocks, as a `"strip"` layout (plots side by side, each
  running the full length) or a `"stack"` layout (two tiers with a
  buffer, halving the width). Treatment order within each block is
  randomised or systematic — a real choice for strip trials, since under
  strong spatial correlation a systematic arrangement can estimate
  contrasts more precisely. The output has the same shape as
  [`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md),
  so sampling locations can be chosen from the plan before any yield
  data exist, and the same code runs before and after harvest.

- The package Title and Description now describe the whole workflow
  rather than the integration step alone. Only four of its functions
  concern integration; the rest design the trial, plan the sampling,
  prepare the data or simulate trials, and the old title made that hard
  to see.

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
