# Subset nucleotide/protein sequences to codon +/- 45 bp context

Subset nucleotide/protein sequences to codon +/- 45 bp context

## Usage

``` r
subset_context(codv, ctx_codons)
```

## Arguments

- codv:

  Annotated variants from \`annotate_coding()\`

- ctx_codons:

  How many flanking codons each to include in the context

## Value

GRanges object with sequence information of only context
