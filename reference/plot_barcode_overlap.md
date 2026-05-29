# Plot barcode overlap between different samples

Plot barcode overlap between different samples

## Usage

``` r
plot_barcode_overlap(all_constructs, valid_barcodes)
```

## Arguments

- all_constructs:

  A named list of all constructs

- valid_barcodes:

  A character vector of possible barcodes (optional)

## Value

A ggplot2 object

## Examples

``` r
if (interactive()) {
    constructs = list(pat1=data.frame(gene_name="GENE1", mut_id="GENE1_A1V",
        pep_id="GENE1_A1V", pep_type="alt", tiled="ATGGCCGCC", barcode_1="AAAA"))
    plot_barcode_overlap(constructs, valid_barcodes=c("AAAA", "CCCC"))
}
```
