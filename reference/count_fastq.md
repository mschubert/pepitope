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
  'sample_id', 'minigenes', 'rep', 'origin', 'barcode' \`minigenes\`
  lists construct libraries separated by \`+\`.

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

A `SummarizedExperiment` object with counts and metadata \`rowData\`
includes \`unassigned_reads\`, the number of reads with a matching
construct barcode but no matching sample barcode.
\`metadata(dset)\$unmapped_reads\` records reads matching neither sample
nor construct barcode, excluding reads too short to extract both
barcodes.

## Examples

``` r
samples = data.frame(sample_id="sample1", minigenes="pat1", rep="1",
    origin="library", barcode="GGG")
constructs = list(pat1=data.frame(gene_name="GENE1", mut_id="GENE1_A1V",
    pep_id="GENE1_A1V", pep_type="alt", tiled="ATGGCCGCC", barcode_1="AAAA"))
fq = tempfile(fileext=".fq")
writeLines(c("@r1", "GGGAAAA", "+", "IIIIIII", "@r2", "GGGAAAA", "+", "IIIIIII"), fq)
count_fastq(fq, samples, constructs, read_structure="3B4M", verbose=FALSE)
#> class: SummarizedExperiment 
#> dim: 1 1 
#> metadata(1): unmapped_reads
#> assays(1): counts
#> rownames(1): AAAA
#> rowData names(8): barcode unassigned_reads ... pep_type tiled
#> colnames(1): sample1
#> colData names(7): sample_id minigenes ... total_reads mapped_reads
```
