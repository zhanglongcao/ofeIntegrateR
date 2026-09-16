# Plot an empirical variogram and the model fitted to it

Points are the empirical semivariances, sized by how many sample pairs
each bin rests on – a bin built from six pairs should not be read as
hard as one built from six hundred, and drawing them the same size is
the commonest way a variogram plot misleads. The fitted model, the
nugget and the practical range are drawn over them.

## Usage

``` r
# S3 method for class 'ofe_variogram'
plot(
  x,
  show_all = FALSE,
  main = NULL,
  xlab = "distance",
  ylab = "semivariance",
  ...
)
```

## Arguments

- x:

  An `ofe_variogram` from
  [`ofe_variogram()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_variogram.md).

- show_all:

  Logical; draw every fitted model rather than only the best.

- main, xlab, ylab:

  Plot labels.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly. Called for the plot.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 30, n_col = 20, n_point_samples = 60,
                          point_range = 6, seed = 1)
pts <- sim$point_samples
pts$x <- pts$col
pts$y <- pts$row
plot(ofe_variogram(pts, value = "point_obs"))

```
