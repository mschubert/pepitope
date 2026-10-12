# Format a read structure from barcode position counts

Raises an error if neither barcode type meets \`min_frac\` or inferred
barcode positions overlap.

## Usage

``` r
.rs_format(counts, sample_width, construct_width, n_reads, min_frac = 0.2)
```

## Arguments

- counts:

  A \`data.frame\` with \`B\`, \`^B\`, \`M\`, and \`^M\` count columns

- sample_width:

  Width of the sample barcode

- construct_width:

  Width of the construct barcode

- n_reads:

  Number of inspected reads

- min_frac:

  Minimum fraction of reads matching each barcode type at its inferred
  position in either orientation

## Value

A read structure string containing the barcode types that meet
\`min_frac\`
