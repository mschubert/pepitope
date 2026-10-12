# Create deterministic example barcodes

Create deterministic example barcodes

## Usage

``` r
example_barcodes(n, width = 12)
```

## Arguments

- n:

  Number of barcodes to create

- width:

  Barcode width in nucleotides

## Value

A character vector of DNA barcodes

## Examples

``` r
example_barcodes(4, width=4)
#> [1] "AAAA" "CAAA" "GAAA" "TAAA"
```
