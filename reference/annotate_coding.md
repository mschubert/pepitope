# Annotate VCF variants with coding changes

Annotate VCF variants with coding changes

## Usage

``` r
annotate_coding(vr, txdb, asm)
```

## Arguments

- vr:

  A VRanges object with SNVs and small indels

- txdb:

  TxDb or EnsDb object

- asm:

  Genomic sequence BSGenome object

## Value

A GRanges object with annotated variants
