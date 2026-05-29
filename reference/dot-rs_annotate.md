# Annotate barcode positions in source FASTQ reads

Annotate barcode positions in source FASTQ reads

## Usage

``` r
.rs_annotate(fq, samples, all_constructs, nrec = 100000L)
```

## Arguments

- fq:

  Path to one FASTQ file

- samples:

  A sample sheet as \`data.frame\` in tsv format. Requires the columns
  'sample_id', 'patient', 'rep', 'origin', 'barcode'

- all_constructs:

  A named list of all construct libraries

- nrec:

  Number of FASTQ records to inspect

## Value

A \`list\` with barcode counts, reads, and inferred read structure
