# Count barcodes using guide-counter

Count barcodes using guide-counter

## Usage

``` r
count_bc(tdir, all_constructs, valid_barcodes, reverse_complement = FALSE)
```

## Arguments

- tdir:

  Path to the directory with demultiplexed FASTQ files

- all_constructs:

  A named list of all construct libraries

- valid_barcodes:

  A character vector of all possible barcodes

- reverse_complement:

  Whether to count the reverse complement of the barcodes instead

## Value

A \`SummarizedExperiment\` object with counts and metadata
