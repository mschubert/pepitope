# Calculate differential abundance of construct barcodes

Calculate differential abundance of construct barcodes

## Usage

``` r
screen_calc(dset, comparisons)
```

## Arguments

- dset:

  A \`SummarizedExperiment\` object from \`count_fastq()\`

- comparisons:

  A character vector of sample and reference condition, or list thereof

## Value

A data.frame of DESeq2 results, or a named list of data.frames

## Examples

``` r
if (interactive()) {
    screen_calc(dset, c("screen", "library"))
}
```
