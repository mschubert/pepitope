# Annotate barcode positions in source FASTQ reads

If structure inference fails, a warning is emitted and the structure is
\`"unknown"\`.

## Usage

``` r
.rs_annotate(fq, samples, all_constructs, nrec = 100000L, min_frac = 0.2)
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

- min_frac:

  Minimum fraction of inspected reads matching each barcode type at its
  inferred position in either orientation

## Value

A \`list\` with barcode counts, reads, and inferred read structure
