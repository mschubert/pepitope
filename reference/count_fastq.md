# Count barcodes directly from source FASTQ files

Count barcodes directly from source FASTQ files

## Usage

``` r
count_fastq(
  fq,
  samples,
  all_constructs,
  valid_barcodes,
  read_structure,
  verbose = TRUE
)
```

## Arguments

- fq:

  Path to one FASTQ file

- samples:

  A sample sheet as `data.frame` in tsv format. Requires the columns
  'sample_id', 'patient', 'rep', 'origin', 'barcode'

- all_constructs:

  A named list of all construct libraries

- valid_barcodes:

  A character vector of all possible construct barcodes

- read_structure:

  A character string describing the FASTQ read structure. If missing,
  this will be inferred from the first reads in \`fq\`.

- verbose:

  Whether to print progress messages (default: TRUE)

## Value

A `SummarizedExperiment` object with counts and metadata

## Examples

``` r
samples = data.frame(sample_id="sample1", patient="pat1", rep="1",
    origin="library", barcode="GGG")
constructs = list(pat1=data.frame(gene_name="GENE1", mut_id="GENE1_A1V",
    pep_id="GENE1_A1V", pep_type="alt", tiled="ATGGCCGCC", barcode_1="AAAA"))
fq = tempfile(fileext=".fq")
writeLines(c("@r1", "GGGAAAA", "+", "IIIIIII", "@r2", "GGGAAAA", "+", "IIIIIII"), fq)
count_fastq(fq, samples, constructs, read_structure="3B4M", verbose=FALSE)
#> class: SummarizedExperiment 
#> dim: 1 1 
#> metadata(0):
#> assays(1): counts
#> rownames(1): AAAA
#> rowData names(7): barcode bc_type ... pep_type tiled
#> colnames(1): sample1
#> colData names(10): sample_id patient ... short label
```
