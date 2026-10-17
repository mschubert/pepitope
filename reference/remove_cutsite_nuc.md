# Remove restriction enzyme sites while preserving amino acids

Uses dynamic programming to minimize nucleotide substitutions while
excluding all recognition sites and their reverse complements everywhere
in the result. Incomplete terminal codons are kept unchanged.
Translation uses the standard genetic code with \`no.init.codon=TRUE\`.

## Usage

``` r
remove_cutsite_nuc(nuc, site, seed = NULL)
```

## Arguments

- nuc:

  A character vector of cDNA nucleotide strings using unambiguous DNA
  bases

- site:

  A character vector of recognition sequences to remove on both strands

- seed:

  Seed for deterministic ordering of equally optimal synonymous codons

## Value

cDNA sequences with the fewest nucleotide substitutions needed to remove
the cut sites while preserving translation
