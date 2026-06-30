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

- **Kriged-covariate** (`fit_integrated_kriged()`): krige the sparse
  point-source variable onto the trial grid with `krige_point_samples()`,
  then include the kriged surface as a fixed covariate in a spatial mixed
  model of the dense response (`ar1(row):ar1(col)` residual via
  [asreml-R](https://vsni.co.uk/software/asreml-r)).
- **Joint bivariate model** (`fit_integrated_joint()`): fit the dense
  response and the point-source measurements *jointly* as two stacked
  response layers sharing a spatial random effect, avoiding a separate
  kriging step.

`simulate_ofe_trial()` generates synthetic test data for both strategies,
and `extract_fixed_effects()` gives a model-agnostic way to pull
treatment-contrast estimates out of either an `lm` or an `asreml` fit.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jeromecy/ofeIntegrateR")
```

The kriged-covariate fallback (`engine = "lm"`) and `krige_point_samples()`
work with no further setup. The full spatial mixed models
(`engine = "asreml"`, and `fit_integrated_joint()`) require
[asreml-R](https://vsni.co.uk/software/asreml-r), a commercial package
from VSNi that is **not** installed automatically — see VSNi for licensing.

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
                                engine = "asreml")  # or engine = "lm"
extract_fixed_effects(fit_a)

# Strategy 2: fit dense + point-source data jointly (no separate kriging step)
fit_b <- fit_integrated_joint(sim$grid, sim$point_samples,
                               response_dense = "dense_response",
                               response_point = "point_obs")
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

## License

MIT © Zhanglong Cao
