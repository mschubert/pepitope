# Use fqtk tool to demultiplex fastq files

Use fqtk tool to demultiplex fastq files

## Usage

``` r
demux_fq(fq, samples, read_structures)
```

## Arguments

- fq:

  A path to the fastq file to demultiplex

- samples:

  A sample sheet as \`data.frame\` in tsv format. Requires the columns
  'sample_id', 'patient', 'rep', 'origin', 'barcode'

- read_structures:

  A character string describing the read structure
