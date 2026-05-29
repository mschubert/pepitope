# Plot annotated read structure examples

Plot annotated read structure examples

## Usage

``` r
plot_read_structure(fq, samples, all_constructs, nrec = 100000L)
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

A \`ggplot2\` object

## Examples

``` r
if (interactive()) {
    samples = data.frame(sample_id="sample1", patient="pat1", rep="1",
        origin="library", barcode="GGG")
    constructs = list(pat1=data.frame(gene_name="GENE1", mut_id="GENE1_A1V",
        pep_id="GENE1_A1V", pep_type="alt", tiled="ATGGCCGCC", barcode_1="AAAA"))
    fq = tempfile(fileext=".fq")
    writeLines(c("@r1", "GGGAAAA", "+", "IIIIIII", "@r2", "GGGAAAA", "+", "IIIIIII"), fq)
    plot_read_structure(fq, samples, constructs, nrec=10)
}
```
