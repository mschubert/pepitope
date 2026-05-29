# Co-culture screen

This guide shows how to calculate differential construct abundance and
visualize the results. To calculate changes in barcode abundance, we use
the `DESeq2` package, commonly used in differential expression analysis.
Finally, we provide a utility function to visualize the results as an
MA-plot (base abundance on the x axis, log2 fold change on the y axis).

## Preparation

In order to perform this analysis, we need the `SummarizedExperiment`
(`dset`) object from the *Quality Control* results, and we need to know
which comparisons we want to perform.

From the *Quality Control* step, we know that we have two repeats each
of `Mock` and `Sample` for a patient-specific (`pat1`) and a common
library (`common`). These two repeats are important, as we need to
estimate within-condition variability as well as between-condition
variability.

However, we need to make sure to perform the screen analysis only on
relevant samples, otherwise the model may mis-estimate the variability:

``` r
dset = dset[,grepl("pat1", dset$patient)]
colData(dset)
#> DataFrame with 4 rows and 10 columns
#>           sample_id     patient       rep      origin     barcode total_reads
#>         <character>    <factor> <numeric> <character> <character>   <numeric>
#> mock1         mock1 pat1+common         1        Mock     TGAGTCC      209496
#> mock2         mock2 pat1+common         2        Mock     CAAGATG      198822
#> screen1     screen1 pat1+common         1      Sample     AACCGAC      407724
#> screen2     screen2 pat1+common         2      Sample     AGAATCG      395252
#>         mapped_reads         smp                short                  label
#>            <numeric> <character>          <character>            <character>
#> mock1         209496      Mock-1   pat1+common Mock-1 pat1+common Mock-1 (..
#> mock2         198822      Mock-2   pat1+common Mock-2 pat1+common Mock-2 (..
#> screen1       407724    Sample-1 pat1+common Sample-1 pat1+common Sample-1..
#> screen2       395252    Sample-2 pat1+common Sample-2 pat1+common Sample-2..
```

## Calculating differential abundance

When we want to perform a comparison between two conditions, we refer in
this comparison to the `origin` column of the sample sheet. In our case,
we have the origin `Mock` and `Sample`, which describe an experiment of
mock-transfected and co-cultured cells and cells transfected with the
actualy construct library, respectively.

Hence, our comparison here is that we want to see the changes of
`Sample` over the `Mock` condition, which we indicate as a character
vector of `c(sample, reference)` or a list thereof.

In case of a single comparison (character vector), a `data.frame` will
be returned. If there are more comparisons supplied in a list of
character vectors, the result will be a list of `data.frame`s:

``` r
res = screen_calc(dset, list(c("Sample", "Mock")))
#> converting counts to integer mode
#> Warning in DESeq2::DESeqDataSet(dset, ~rep + origin): some variables in design
#> formula are characters, converting to factors
#>   the design formula contains one or more numeric variables with integer values,
#>   specifying a model with increasing fold change for higher values.
#>   did you mean for this to be a factor? if so, first convert
#>   this variable to a factor using the factor() function
#> using pre-existing size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> -- note: fitType='parametric', but the dispersion trend was not well captured by the
#>    function: y = a/x + b, and a local regression fit was automatically substituted.
#>    specify fitType='local' or 'mean' to avoid this message next time.
#> final dispersion estimates
#> fitting model and testing
#> Joining with `by = join_by(barcode)`
```

The result will look like the following:

``` r
res[[1]]
#> # A tibble: 206 × 19
#>    barcode baseMean log2FoldChange lfcSE   stat   pvalue     padj bc_type var_id
#>    <chr>      <dbl>          <dbl> <dbl>  <dbl>    <dbl>    <dbl> <fct>   <chr> 
#>  1 AACAAC…    1740.         -3.96  0.199 -19.9  7.10e-88 1.46e-85 pat1    chr1:…
#>  2 AACGCC…    3450.          0.662 0.192   3.44 5.71e- 4 5.89e- 2 common  chr18…
#>  3 AACACA…    1104.          0.502 0.167   3.00 2.68e- 3 1.84e- 1 pat1    chr15…
#>  4 AACAAG…    2306.          0.626 0.220   2.84 4.52e- 3 1.84e- 1 pat1    chr7:…
#>  5 AACACG…    1291.         -0.438 0.156  -2.82 4.85e- 3 1.84e- 1 pat1    chr4:…
#>  6 AACGAT…    1042.          0.488 0.176   2.78 5.47e- 3 1.84e- 1 common  chr1:…
#>  7 AACAAG…    2660.          0.423 0.158   2.68 7.42e- 3 1.84e- 1 pat1    chr7:…
#>  8 AACGAC…    2155.          0.426 0.160   2.66 7.87e- 3 1.84e- 1 common  chr9:…
#>  9 AACAGA…    1852.         -0.461 0.174  -2.65 8.12e- 3 1.84e- 1 pat1    chr9:…
#> 10 AACGAT…    2303.         -0.421 0.162  -2.60 9.38e- 3 1.84e- 1 common  ETV6-…
#> # ℹ 196 more rows
#> # ℹ 10 more variables: mut_id <chr>, pep_id <chr>, pep_type <chr>,
#> #   gene_name <chr>, gene_id <chr>, tx_id <chr>, tiled <chr>, n_tiles <dbl>,
#> #   nt <dbl>, peptide <chr>
```

## Plotting the screen results

We can plot the result of this differential abundance analysis using the
`plot_screen` function. The only argument we need supply is a results
table, but we can fine-tune the plot with the following additional
arguments:

- `sample` – which library (or libraries) to plot (default: all)
- `links` – whether to draw arrows between ref and significant alt
  peptides
- `labs` – whether to label constructs in less dense areas of the plot
- `cap_fc` – a maximum (and negative minimum) fold change value to bound
  data points to

``` r
plot_screen(res$`Sample vs Mock`)
#> Joining with `by = join_by(bc_type, mut_id)`
#> Registered S3 methods overwritten by 'ggpp': method from
#> heightDetails.titleGrob ggplot2 widthDetails.titleGrob ggplot2
```

![](screen_files/figure-html/unnamed-chunk-7-1.png)

As with the *Quality Control* plots already, we can also create and
interactive plot. This will not display links between `ref` and `alt`
peptides, but instead highlight a contruct group (by `mut_id`) when
interacting with a data point.
