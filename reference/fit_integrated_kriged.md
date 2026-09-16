# Fit a treatment model using a kriged point-source covariate

Implements the "kriged-covariate" data integration strategy: a sparse
point-source variable (already kriged onto the trial grid, e.g. via
[`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md))
is included as a fixed covariate alongside the treatment factor, with an
AR1xAR1 spatial residual structure to absorb any remaining spatial
autocorrelation in the dense response.

## Usage

``` r
fit_integrated_kriged(
  data,
  response,
  treat,
  covariate = NULL,
  random = NULL,
  row = "row",
  col = "col",
  engine = c("ofe", "asreml", "lme", "gls", "lm"),
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

- random:

  Optional one-sided formula of random effects, such as `~ rep` or
  `~ block`. Real strip trials are replicated, so this is usually
  needed. Supported by `engine = "ofe"`, `"asreml"` and `"lme"`; the
  `"gls"` and `"lm"` engines cannot fit random effects and raise an
  error rather than ignoring the argument.

- row, col:

  Character; names of the row/column position columns used to build the
  spatial residual structure. Ignored when `engine = "lm"`. For a
  pseudo-environment residual, or any structure other than
  `ar1(row):ar1(col)`, call
  [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
  directly with your own `residual` formula.

- engine:

  Character; `"ofe"` (default) fits `response ~ treat + covariate` with
  an `ar1(row):ar1(col)` residual by REML via
  [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
  — the same model asreml fits, written in base R and needing no
  licence. `"asreml"` fits the identical model through
  `asreml::asreml()` where a licence is available; the two agree to
  several significant figures on the structures both support, and asreml
  is faster on large lattices. `"lme"` fits random effects **and** an
  exponential spatial correlation via
  [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html), which is the
  open-source counterpart to the asreml fit for a replicated trial.
  `"gls"` fits the same fixed-effects model via
  [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html) with an
  exponential spatial correlation structure
  ([`nlme::corExp()`](https://rdrr.io/pkg/nlme/man/corExp.html) on
  `row`/`col`) — an open-source (CRAN-only, no licence required)
  alternative that still accounts for residual spatial autocorrelation.
  `"lm"` fits the same fixed-effects model via
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html) with no spatial
  residual structure at all — the simplest fallback, but does not
  account for residual spatial autocorrelation.

- ...:

  Additional arguments passed to the engine:
  [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
  (e.g. `control = ofe_control(trace = TRUE)`), `asreml::asreml()` (e.g.
  `maxit`), [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html),
  [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html), or
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html).

## Value

The fitted model object (class `ofe_fit`, `asreml`, `lme`, `gls`, or
`lm`). Use
[`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
to retrieve a tidy table of fixed-effect estimates, including the
treatment contrasts.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 12, seed = 1)
krieged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
#> Warning: No convergence after 200 iterations: try different initial values?
fit <- fit_integrated_kriged(krieged, response = "dense_response",
                              treat = "treat", covariate = "point_obs_kriged",
                              row = "row", col = "col")
extract_fixed_effects(fit)
#>               term    estimate        se
#> 1      (Intercept) -0.09683842 0.2851930
#> 2           treatB  1.01542575 0.1486771
#> 3           treatC  1.85459880 0.1453960
#> 4 point_obs_kriged  0.53862942 0.2013503

if (requireNamespace("asreml", quietly = TRUE)) {
  fit_asr <- fit_integrated_kriged(krieged, response = "dense_response",
                                    treat = "treat", covariate = "point_obs_kriged",
                                    row = "row", col = "col", engine = "asreml")
  extract_fixed_effects(fit_asr)
}
```
