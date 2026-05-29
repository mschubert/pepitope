# Remove a Restriction Enzyme cut site but keep AA in a tiled peptide data.frame

Remove a Restriction Enzyme cut site but keep AA in a tiled peptide
data.frame

## Usage

``` r
remove_cutsite(pep, ...)
```

## Arguments

- pep:

  A data.frame of tiled peptides

- ...:

  Named argumennts of cut sites, e.g. \`BbsI="GAAGAC"\`

## Value

A data.frame with replace nucleotides and number of replacements

## Examples

``` r
pep = data.frame(pep_id="p1", tiled="ATGGCCGCC")
remove_cutsite(pep, BbsI="GAAGAC")
#>   pep_id     tiled BbsI_replaced
#> 1     p1 ATGGCCGCC             0
```
