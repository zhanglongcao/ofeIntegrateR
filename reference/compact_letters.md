# Compact letter display from any set of pairwise comparisons

Turns pairwise p-values into the a/b/c annotation: two levels share a
letter when they are not significantly different. Levels are ordered by
their means where these are supplied, so `"a"` marks the best treatment.

## Usage

``` r
compact_letters(
  comparisons,
  alpha = 0.05,
  means = NULL,
  decreasing = TRUE,
  p_col = "p.value"
)
```

## Arguments

- comparisons:

  Either a symmetric matrix of p-values whose dimnames are the level
  names, or a data frame with a p-value column and the two levels of
  each comparison. The level columns are taken from `level1` and
  `level2` if present, otherwise from a `contrast` column of the form
  `"B - A"`, otherwise from the first two non-numeric columns.

- alpha:

  Significance level; comparisons with `p.value > alpha` are treated as
  "not different" and so share a letter.

- means:

  Optional named numeric vector of the means, used to order the levels.
  Without it the levels keep the order they appear in.

- decreasing:

  Logical; with `means` supplied, order from the largest mean down, so
  `"a"` goes to the highest. Set `FALSE` if smaller is better.

- p_col:

  Character; name of the p-value column when `comparisons` is a data
  frame.

## Value

A named character vector of letter groups, one per level, in the order
the letters were assigned.

## Details

This is deliberately separate from
[`ofe_lsd()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_lsd.md)
and takes only the comparisons, so it can letter the output of anything
– `asreml::predict()`, `emmeans`, a table typed in by hand – and not
only a fit from this package.

## See also

[`ofe_lsd()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_lsd.md)
to go straight from a fitted model to a lettered table.

## Examples

``` r
cmp <- data.frame(
  level1 = c("A", "A", "B"),
  level2 = c("B", "C", "C"),
  p.value = c(0.40, 0.001, 0.002)
)
compact_letters(cmp, means = c(A = 3.0, B = 3.1, C = 5.0))
#>   C   B   A 
#> "a" "b" "b" 
```
