# Fit a treatment model using a kriged point-source covariate

Implements the "kriged-covariate" data integration strategy: a sparse
point-source variable (already kriged onto the trial grid, e.g. via
[`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md))
is included as a fixed covariate alongside the treatment factor, with an
AR1xAR1 spatial residual structure (when `engine = "asreml"`) to absorb
any remaining spatial autocorrelation in the dense response.

## Usage

``` r
fit_integrated_kriged(
  data,
  response,
  treat,
  covariate = NULL,
  row = "row",
  col = "col",
  engine = c("asreml", "gls", "lm"),
  ...
)
```

## Arguments

- data:

  Data frame containing the response, treatment factor, kriged
  covariate, and row/column position columns (one row per grid cell).

- response:

  Character; name of the dense response column (e.g. yield).

- treat:

  Character; name of the treatment factor column.

- covariate:

  Character vector; name(s) of the kriged point-source covariate
  column(s) (e.g. the `*_kriged` columns produced by
  [`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md)).
  Supply several names to integrate more than one point variable at
  once, such as soil nitrogen and phosphorus. Use `NULL` to fit the
  treatment model with no point-source covariate at all: this is the
  baseline that an integrated analysis has to beat, since it uses only
  the dense layer and costs no sampling.

- row, col:

  Character; names of the row/column position columns used to build the
  spatial residual structure. Ignored when `engine = "lm"`.

- engine:

  Character; `"asreml"` (default) fits `response ~ treat + covariate`
  with an `ar1(row):ar1(col)` residual via `asreml::asreml()` — requires
  a licensed copy of asreml-R. `"gls"` fits the same fixed-effects model
  via [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html) with an
  exponential spatial correlation structure
  ([`nlme::corExp()`](https://rdrr.io/pkg/nlme/man/corExp.html) on
  `row`/`col`) — an open-source (CRAN-only, no licence required)
  alternative that still accounts for residual spatial autocorrelation.
  `"lm"` fits the same fixed-effects model via
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html) with no spatial
  residual structure at all — the simplest fallback, but does not
  account for residual spatial autocorrelation.

- ...:

  Additional arguments passed to `asreml::asreml()` (e.g. `maxit`),
  [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html), or
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html).

## Value

The fitted model object (class `asreml`, `gls`, or `lm`). Use
[`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
to retrieve a tidy table of fixed-effect estimates, including the
treatment contrasts.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 12, seed = 1)
krieged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
#> Warning: No convergence after 200 iterations: try different initial values?
fit_lm <- fit_integrated_kriged(krieged, response = "dense_response",
                                 treat = "treat", covariate = "point_obs_kriged",
                                 row = "row", col = "col", engine = "lm")
extract_fixed_effects(fit_lm)
#>               term    estimate        se
#> 1      (Intercept) -0.07102343 0.2878442
#> 2           treatB  1.00366194 0.1133115
#> 3           treatC  1.90258169 0.1060346
#> 4 point_obs_kriged  0.57143055 0.2116354

if (requireNamespace("asreml", quietly = TRUE)) {
  fit_asr <- fit_integrated_kriged(krieged, response = "dense_response",
                                    treat = "treat", covariate = "point_obs_kriged",
                                    row = "row", col = "col", engine = "asreml")
  extract_fixed_effects(fit_asr)
}
```
