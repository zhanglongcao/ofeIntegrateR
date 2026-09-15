# Fit a joint bivariate spatial model for dense and point-source data

Implements the "joint model" data integration strategy: rather than
pre-kriging the point-source variable and using it as a fixed covariate
(see
[`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)),
the dense response and the point-source measurements are fitted
**jointly** as two stacked response "layers" sharing a spatial random
effect at each grid location. This lets sparse point samples directly
inform the dense response's spatial structure (and vice versa) without a
separate kriging step.

## Usage

``` r
fit_integrated_joint(
  grid,
  point_samples,
  response_dense,
  response_point,
  treat = "treat",
  row = "row",
  col = "col",
  engine = c("asreml", "sommer"),
  cor_structure = c("independent", "unstructured"),
  min_point_n_for_us = 30,
  maxit = 60,
  ai_sing = TRUE,
  ...
)
```

## Arguments

- grid:

  Data frame, one row per grid location, with the dense response,
  treatment factor, and row/column position columns (as produced by
  [`simulate_ofe_trial()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_ofe_trial.md)'s
  `grid` element).

- point_samples:

  Data frame of sparse point-source observations with row/column
  position columns and the point-source value column (as produced by
  [`simulate_ofe_trial()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_ofe_trial.md)'s
  `point_samples` element).

- response_dense:

  Character; name of the dense response column in `grid` (e.g. yield).

- response_point:

  Character; name of the point-source value column in `point_samples`
  (e.g. soil, tissue, or disease measurement).

- treat:

  Character; name of the treatment factor column in `grid`. The
  treatment effect is estimated from the dense layer only.

- row, col:

  Character; names of the row/column position columns, shared between
  `grid` and `point_samples`.

- engine:

  Character; `"asreml"` (default) fits the model via `asreml::asreml()`
  — requires a licensed copy of asreml-R. `"sommer"` fits an open-source
  (CRAN-only, no licence required) analogue via
  [`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html): the
  dense and point-source responses are modelled as two traits in a
  multi-trait model sharing a random effect per grid location
  (`vsr(unit, Gtc = ...)`), with the treatment effect estimated
  separately for each trait (only the dense-trait estimate is
  meaningful; see Details).

- cor_structure:

  Character; `"independent"` (default) fits `diag(layer):id(unit)`
  (asreml) or `vsr(unit, Gtc = diag(2))` (sommer), i.e. separate spatial
  variances for the dense and point-source layers with no
  cross-covariance — robust when point samples are sparse.
  `"unstructured"` fits `us(layer):id(unit)` (asreml) or
  `vsr(unit, Gtc = unsm(2))` (sommer), which additionally estimates the
  spatial cross-covariance between layers (how strongly point-source
  spatial pattern informs the dense response), but is often
  unidentifiable below `min_point_n_for_us` observations.

- min_point_n_for_us:

  Integer; minimum number of non-missing point-source observations below
  which a warning is issued when `cor_structure = "unstructured"`.

- maxit:

  Maximum number of asreml iterations. Ignored when `engine = "sommer"`.

- ai_sing:

  Logical; if `TRUE` (default), temporarily set the session-level
  `asreml.options(ai.sing = TRUE)` for the duration of this fit
  (restored to its previous value on exit), so that singularities in the
  Average Information matrix do not abort estimation. In testing, the
  layer-sharing random effect (even with
  `cor_structure = "independent"`) can trigger this with sparse
  point-source data, so it is enabled by default. Ignored when
  `engine = "sommer"`.

- ...:

  Additional arguments passed to `asreml::asreml()` or
  [`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html).

## Value

A fitted `asreml` or `mmer` model object. Treatment contrasts for the
dense layer can be retrieved with
[`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
(look for terms matching `at(layer, "dense"):<treat>` for asreml, or
`<response_dense>:<treat><level>` for sommer).

## Details

For `engine = "asreml"`, residuals are modelled with
`dsum(~units | layer)`, i.e. separate i.i.d. residual variances per
layer. An `ar1(row):ar1(col)` residual was found to be unidentifiable
for the point-source layer in testing, because a handful of sparse,
irregularly placed observations cannot support a full spatial
autocorrelation structure on their own; the shared random effect
(`diag(layer)` or `us(layer)`) is what carries spatial information
between layers instead.

For `engine = "sommer"`, the same idea is expressed as a multi-trait
model (`cbind(dense, point) ~ treat`) with a shared random effect per
grid location (`unit`) and heterogeneous, uncorrelated residual
variances per trait (`vsr(units, Gtc = diag(2))`). sommer estimates the
treatment effect separately for *each* trait by default; the
point-source trait's treatment estimate is not meaningful (point samples
carry no treatment information of their own) and should be ignored —
only the dense-trait estimate is the quantity of interest. Rows are kept
even when only one of the two responses is observed via
`naMethodY = "include"`. In testing against simulated data, this
sommer-based fit reproduced the asreml `cor_structure = "independent"`
treatment-contrast estimates almost exactly.

## Examples

``` r
if (requireNamespace("asreml", quietly = TRUE)) {
  sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 12, seed = 1)
  fit <- fit_integrated_joint(sim$grid, sim$point_samples,
                               response_dense = "dense_response",
                               response_point = "point_obs")
  extract_fixed_effects(fit)
}

if (requireNamespace("sommer", quietly = TRUE)) {
  # sommer's multi-trait model needs more point-source observations than
  # asreml's to avoid a singular system; a sparser grid like the asreml
  # example above is not guaranteed to converge.
  sim <- simulate_ofe_trial(n_row = 30, n_col = 21, n_point_samples = 25, seed = 11)
  fit <- fit_integrated_joint(sim$grid, sim$point_samples,
                               response_dense = "dense_response",
                               response_point = "point_obs",
                               engine = "sommer")
  extract_fixed_effects(fit)
}
#>                         term    estimate         se
#> 1 dense_response:(Intercept)  0.39432617 0.04629048
#> 2      point_obs:(Intercept)  0.14632979 0.01408539
#> 3      dense_response:treatB  0.15606508 0.06546463
#> 4           point_obs:treatB -0.04129616 0.01991974
#> 5      dense_response:treatC  1.31882933 0.06546463
#> 6           point_obs:treatC -0.02382778 0.01991974
```
