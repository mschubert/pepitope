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

## Examples

``` r
if (interactive()) {
    res = data.frame(baseMean=c(20, 30), log2FoldChange=c(-1, 1), padj=c(0.2, 0.05),
        gene_name=c("GENE1", "GENE1"), pep_type=c("ref", "alt"),
        bc_type=c("pat1", "pat1"), mut_id=c("GENE1_A1", "GENE1_A1"))
    plot_screen(res, links=FALSE, labs=FALSE)
}
```
