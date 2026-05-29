# Subset the peptide context for gene fusions

Subset the peptide context for gene fusions

## Usage

``` r
subset_context_fusion(res, ctx_codons)
```

## Arguments

- res:

  A DataFrame object from \`fusions\`

- ctx_codons:

  How many flanking codons each to include in the context

## Value

A DataFrame object of gene fusions

## Examples

``` r
if (interactive()) {
    subset_context_fusion(fusions, ctx_codons=15)
}
```
