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

  Seed for deterministic pseudo-random candidate ordering

## Value

cDNA with minimal changes to no longer contain the cut site
