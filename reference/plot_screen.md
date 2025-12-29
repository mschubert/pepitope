# Plot screen results

Plot screen results

## Usage

``` r
plot_screen(res, sample = NULL, links = TRUE, labs = TRUE, cap_fc = 8)
```

## Arguments

- res:

  A results \`data.frame\` from \`screen_calc()\`

- sample:

  Which library to plot (default: all)

- links:

  Whether to draw arrows between ref and significant alt peptides

- labs:

  Whether to label genes in less dense areas

- cap_fc:

  Maximum amount of fold-change to limit values to

## Value

A \`ggplot2\` object of the differential expression results
