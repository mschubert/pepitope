# Use \`guide-counter\` via a system call to actually count

Use \`guide-counter\` via a system call to actually count

## Usage

``` r
count_external(
  tdir,
  valid_barcodes,
  reverse_complement = FALSE,
  exact_match = TRUE
)
```

## Arguments

- tdir:

  Path to the directory with demultiplexed FASTQ files

- valid_barcodes:

  A character vector of all possible barcodes

- reverse_complement:

  Whether to count the reverse complement of the barcodes instead
  (default: FALSE)

- exact_match:

  Whether the read needs to contain the exact barcode (default: TRUE)

## Value

A list with the data.frame meta and matrix counts
