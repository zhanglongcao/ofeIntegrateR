# ofeIntegrateR

<!-- badges: start -->
[![R-CMD-check](https://github.com/jeromecy/ofeIntegrateR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jeromecy/ofeIntegrateR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Tools for integrating sparse point-source measurements (e.g. soil cores,
tissue samples, disease ratings) with dense spatial covariates (e.g. EM38
surveys, yield maps) in on-farm experimentation (OFE) strip trials.

Developed for Milestone 4 of the AAGI-CU-RD-OFE GRDC project ("Development
of processes to integrate point-source data and high-resolution data").

## Two integration strategies

Both strategies have a fully open-source path (no licence required) as
well as an [asreml-R](https://vsni.co.uk/software/asreml-r) path for users
with a licence:

- **Kriged-covariate** (`fit_integrated_kriged()`): krige the sparse
  point-source variable onto the trial grid with `krige_point_samples()`,
  then include the kriged surface as a fixed covariate in a spatial model
  of the dense response. `engine = "asreml"` fits an `ar1(row):ar1(col)`
  residual via asreml-R; `engine = "gls"` (open source) fits an
  exponential spatial correlation structure via `nlme::gls()`;
  `engine = "lm"` (open source) ignores spatial autocorrelation entirely
  — the simplest but weakest fallback.
- **Joint bivariate model** (`fit_integrated_joint()`): fit the dense
  response and the point-source measurements *jointly* as two response
  layers/traits sharing a spatial random effect, avoiding a separate
  kriging step. `engine = "asreml"` fits this via asreml-R;
  `engine = "sommer"` (open source) fits an analogous multi-trait model
  via `sommer::mmer()`. In testing, the sommer fit reproduced the asreml
  treatment-contrast estimates almost exactly.

`simulate_ofe_trial()` generates synthetic test data for both strategies,
and `extract_fixed_effects()` gives a model-agnostic way to pull
treatment-contrast estimates out of an `lm`, `gls`, `asreml`, or `mmer` fit.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jeromecy/ofeIntegrateR")
```

Everything works out of the box with open-source dependencies only
(`engine = "lm"` / `"gls"` for `fit_integrated_kriged()`, and
`engine = "sommer"` for `fit_integrated_joint()`). The `engine = "asreml"`
options additionally require [asreml-R](https://vsni.co.uk/software/asreml-r),
a commercial package from VSNi that is **not** installed automatically —
see VSNi for licensing.

## Example

```r
library(ofeIntegrateR)

# Simulate a strip trial with a sparse point-source covariate
sim <- simulate_ofe_trial(n_row = 30, n_col = 21, n_treat = 3,
                           treat_effects = c(0, 0.8, 1.6),
                           n_point_samples = 25, seed = 11)

# Strategy 1: krige the point-source variable, use it as a fixed covariate
kriged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
fit_a <- fit_integrated_kriged(kriged, response = "dense_response",
                                treat = "treat", covariate = "point_obs_kriged",
                                engine = "gls")  # or "asreml" / "lm"
extract_fixed_effects(fit_a)

# Strategy 2: fit dense + point-source data jointly (no separate kriging step)
fit_b <- fit_integrated_joint(sim$grid, sim$point_samples,
                               response_dense = "dense_response",
                               response_point = "point_obs",
                               engine = "sommer")  # or "asreml"
extract_fixed_effects(fit_b)
```

## Notes on asreml-R behaviour

A few asreml-R quirks this package works around, documented here so they
aren't rediscovered the hard way:

- `summary(fit)$coef.fixed` returns `NULL` in asreml-R >= 4.2; you need
  `summary(fit, coef = TRUE)$coef.fixed`. `extract_fixed_effects()`
  handles this for you.
- `ai.sing` is a *session-level* option (`asreml::asreml.options(ai.sing = TRUE)`),
  not an argument to `asreml()` — passing `ai.sing = TRUE` directly to
  `asreml()` is silently ignored. `fit_integrated_joint()` sets it
  temporarily (restored on exit) via the `ai_sing` argument (default
  `TRUE`), since the shared-layer random effect can trigger Average
  Information matrix singularities even with sparse, well-behaved data.
- An `ar1(row):ar1(col)` residual is generally unidentifiable for the
  point-source layer in the joint model, because a handful of sparse,
  irregularly placed observations cannot support a full spatial
  autocorrelation structure on their own. `fit_integrated_joint()` instead
  uses `dsum(~units | layer)` (separate i.i.d. residual variances per
  layer) and lets the shared random effect (`diag(layer)` or `us(layer)`)
  carry spatial information between layers.

## Notes on the sommer fallback

- `sommer::mmer()` defaults to `naMethodY = "exclude"`, which drops a row
  if *either* trait is missing — catastrophic here, since the
  point-source trait is `NA` almost everywhere. `fit_integrated_joint()`
  sets `naMethodY = "include"` so the dense response keeps using all its
  rows.
- sommer fits the treatment effect separately for both the dense and
  point-source traits; only the dense-trait estimate is meaningful and
  should be used (`extract_fixed_effects()` returns terms for both, named
  e.g. `dense_response:treatB` and `point_obs:treatB` — ignore the latter).
- The sommer multi-trait model needs enough point-source observations
  relative to the grid size to avoid a singular system (`vsr(unit, ...)`
  fits one random-effect level per grid cell). A sparse setup that works
  fine in asreml (e.g. 200 grid cells, 12 points) can fail outright in
  sommer; a denser one (e.g. 630 grid cells, 25 points) works reliably in
  testing.

## License

MIT © Zhanglong Cao
