# Remove restriction enzyme sites while preserving amino acids

Minimizes nucleotide substitutions while preserving translation and
excluding all supplied recognition sites on both strands, including
sites created by other edits.

## Usage

``` r
remove_cutsite(pep, ...)
```

## Arguments

- pep:

  A data.frame of tiled peptides

- ...:

  Named recognition sequences containing A, C, G, and T, e.g.
  \`BbsI="GAAGAC", EcoRI="GAATTC"\`. All sites are removed together.

## Value

A data.frame with synonymously recoded sequences and the number of
original occurrences on either strand in one \`\<name\>\_replaced\`
column per site. Palindromic recognition sites are counted once per
occurrence.

## Examples

``` r
pep = data.frame(pep_id="p1", tiled="ATGGAAGACGAATTC")
remove_cutsite(pep, BbsI="GAAGAC", EcoRI="GAATTC")
#>   pep_id           tiled BbsI_replaced EcoRI_replaced
#> 1     p1 ATGGAGGACGAGTTC             1              1
```
