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

  Seed for deterministic pseudo-random count generation

## Value

The path to the created FASTQ file

## Examples

``` r
samples = data.frame(sample_id="sample1", patient="pat1", rep="1",
    origin="library", barcode="GGG")
constructs = list(pat1=data.frame(gene_name="GENE1", mut_id="GENE1_A1V",
    pep_id="GENE1_A1V", pep_type="alt", tiled="ATGGCCGCC", barcode_1="AAAA"))
example_fastq(samples, constructs, target_reads=2, custom=FALSE)
#> [1] "/tmp/RtmphuvFfx/my_seqdata.fq"
```
