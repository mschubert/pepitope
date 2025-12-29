# Get the alternative coding sequence of a variant or gene fusion

This takes into account any (newly introduced) STOP codons and UTR
readthroughs

## Usage

``` r
get_coding_seq(asm, txdb, ..., include_stop = TRUE)
```

## Arguments

- asm:

  Genomic sequence BSGenome object

- txdb:

  A transcription database, eg. AnnotationHub()\[\["AH100643"\]\]

- ...:

  A named DNAStringSet objects where each row is translated
  consecutively

- include_stop:

  Whether to include the STOP codon

## Value

A merged DNAStringSet object with the translated nucleotides
