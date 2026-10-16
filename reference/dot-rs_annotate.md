# Annotate barcode positions in source FASTQ reads

If structure inference fails, a warning is emitted and the structure is
\`"unknown"\`.

## Usage

``` r
.rs_annotate(
  fq,
  sample_barcodes,
  construct_barcodes,
  nrec = 100000L,
  min_frac = 0.2
)
```

## Arguments

- fq:

  Path to one FASTQ file

- sample_barcodes:

  A character vector of known sample barcodes

- construct_barcodes:

  A character vector of known construct barcodes

- nrec:

  Number of FASTQ records to inspect

- min_frac:

  Minimum fraction of inspected reads matching each barcode type at its
  inferred position in either orientation

## Value

A \`list\` with barcode counts, reads, and inferred read structure
