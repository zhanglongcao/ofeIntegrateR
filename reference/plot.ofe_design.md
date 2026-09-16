# Draw the trial layout

A plan of the plots, coloured and labelled by treatment. Look at it
before the trial goes in: a randomisation that happens to put the same
rate at both ends of a slope, or a buffer in the wrong place, is obvious
on the map and invisible in the design table.

## Usage

``` r
# S3 method for class 'ofe_design'
plot(x, label = TRUE, rep_borders = TRUE, main = NULL, ...)
```

## Arguments

- x:

  An `ofe_design` from
  [`make_trial_design()`](https://www.zcao.space/ofeIntegrateR/reference/make_trial_design.md).

- label:

  Logical; write the treatment on each plot.

- rep_borders:

  Logical; outline the replicate blocks.

- main:

  Plot title.

- ...:

  Passed to
  [`ofe_map()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_map.md).

## Value

`x`, invisibly.

## Details

Every plot carries its treatment label as well as its colour, so the
plan stays readable in greyscale, to a colour-blind reader, and past the
six treatments the categorical palette is validated for.

## Examples

``` r
d <- make_trial_design(c("N0", "N60", "N120"), n_rep = 4, seed = 1)
plot(d)

```
