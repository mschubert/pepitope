# Simulate sequencing data and write them to a temporary FASTQ file

Simulate sequencing data and write them to a temporary FASTQ file

## Usage

``` r
example_fastq(
  samples,
  peptide_sheets,
  target_reads = 1000,
  custom = TRUE,
  seed = 91651
)
```

## Arguments

- samples:

  The .tsv or data.frame file containing sample information

- peptide_sheets:

  A list, each item containing construct information

- target_reads:

  How many reads to simulate on average

- custom:

  Whether to add custom modifications to founds

- seed:

  The random seed used for sampling the number of reads

## Value

The path to the created FASTQ file
