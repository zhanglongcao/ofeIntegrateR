# Predicted means with LSD letters

The table a trial report prints: each level's predicted mean, its
standard error, and the a/b/c letters saying which means are separable.
Levels sharing a letter are not significantly different at `alpha`.

## Usage

``` r
ofe_lsd(
  object,
  term,
  by = NULL,
  alpha = 0.05,
  adjust = "none",
  use = c("pairwise", "lsd"),
  decreasing = TRUE,
  sort = TRUE,
  level = 0.95
)
```

## Arguments

- object:

  An `ofe_fit` from
  [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md).

- term:

  Character; the factor whose means are to be lettered.

- by:

  Character or `NULL`; letter the means within each level of a second
  factor, as in `ofe_lsd(fit, "treat", by = "zone")` for a
  pseudo-environment analysis. Comparisons and letters stay inside a
  level; the same letter in two zones means nothing.

- alpha:

  Significance level for the letters (default 0.05).

- adjust:

  Multiplicity adjustment for the pairwise p-values: `"none"` (default),
  `"tukey"`, `"sidak"`, or any method taken by
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Ignored
  when `use = "lsd"`.

- use:

  `"pairwise"` to judge each comparison on its own SED, or `"lsd"` to
  apply a single average LSD to all of them (see above).

- decreasing:

  Logical; assign `"a"` to the largest mean (default). Set `FALSE` where
  smaller is better, such as a disease score.

- sort:

  Logical; order the rows by mean rather than by factor level.

- level:

  Numeric; confidence level for the interval on each mean.

## Value

A data frame with the `by` level where given, the `term` level,
`estimate`, `std.error`, `lower`, `upper` and `group` (the letters).
Attached as attributes: `lsd`, a one-row-per-group data frame of the
average SED, the LSD, the degrees of freedom and the settings used; and
`comparisons`, the full pairwise table behind the letters.

## Details

Two ways of deciding "different" are offered, and they part company on
exactly the data this package is for.

- `use = "pairwise"` (default):

  Each pair is judged on its own standard error of difference. This is
  the right choice for a spatial model of a strip trial: neighbouring
  strips are compared more precisely than distant ones, and the SEDs
  genuinely differ.

- `use = "lsd"`:

  One least significant difference, \\t\_{1-\alpha/2, \nu} \times
  \overline{SED}\\, applied to every comparison. This is the textbook
  LSD and what a balanced design implies; it is also what most published
  tables mean by "LSD (5%)". Quote the value with the table.

The default `adjust = "none"` gives unprotected LSD comparisons, which
is the convention these letters normally carry – and which does not
control the family-wise error rate. With more than three or four
treatments, say so in the caption or pass `adjust = "tukey"`.

## See also

[`ofe_means()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_means.md)
for the means and comparisons on their own, and
[`compact_letters()`](https://www.zcao.space/ofeIntegrateR/reference/compact_letters.md)
to letter comparisons from any other source.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 20, n_col = 12, seed = 1)
fit <- fit_ofe(dense_response ~ treat, data = sim$grid)

tab <- ofe_lsd(fit, "treat")
tab
#>   treat  estimate std.error     lower    upper group
#> 1     C 1.6736653 0.1251896 1.4270388 1.920292     a
#> 2     B 1.5182691 0.1245343 1.2729335 1.763605     a
#> 3     A 0.8699961 0.1251896 0.6233696 1.116623     b
attr(tab, "lsd")
#>   average_sed      lsd  df alpha adjust      use
#> 1   0.1720583 0.338959 237  0.05   none pairwise

# Protected against multiplicity instead
ofe_lsd(fit, "treat", adjust = "tukey")
#>   treat  estimate std.error     lower    upper group
#> 1     C 1.6736653 0.1251896 1.4270388 1.920292     a
#> 2     B 1.5182691 0.1245343 1.2729335 1.763605     a
#> 3     A 0.8699961 0.1251896 0.6233696 1.116623     b
```
