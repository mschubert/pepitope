# Create example peptide sheets for multiple samples

Create example peptide sheets for multiple samples

## Usage

``` r
example_peptides(valid_barcodes)
```

## Arguments

- valid_barcodes:

  A character vector of valid barcodes

## Value

A named list of peptide/minigene constructs with barcodes

## Examples

``` r
bases = expand.grid(rep(list(c("A", "C", "G", "T")), 4))
valid_barcodes = apply(bases, 1, paste0, collapse="")
if (interactive()) {
    example_peptides(valid_barcodes)
}
```
