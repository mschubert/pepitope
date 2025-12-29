# Remove a Restriction Enzyme cut site but keep AA

Remove a Restriction Enzyme cut site but keep AA

## Usage

``` r
remove_cutsite_nuc(nuc, site, seed = NULL)
```

## Arguments

- nuc:

  cDNA nucleotide string

- site:

  Recognition site to be replaced (fwd+rev comp)

- seed:

  Set random seed to select same changes on multiple runs

## Value

cDNA with minimal changes to no longer contain the cut site
