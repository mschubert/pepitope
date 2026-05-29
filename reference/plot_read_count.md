# Plot the overall read counts

Plot the overall read counts

## Usage

``` r
plot_read_count(dset, n_cutoff = 10)

plot_reads(dset)
```

## Arguments

- dset:

  The \`SummarizedExperiment\` object from \`count_fastq\`

- n_cutoff:

  Minimum reads for counting a barcode as represented

## Value

A patchwork object containing the read count summary plots

## Examples

``` r
if (interactive()) {
    samples = data.frame(sample_id="sample1", patient="pat1", rep="1",
        origin="library", barcode="GGG")
    constructs = list(pat1=data.frame(gene_name="GENE1", mut_id="GENE1_A1V",
        pep_id="GENE1_A1V", pep_type="alt", tiled="ATGGCCGCC", barcode_1="AAAA"))
    fq = tempfile(fileext=".fq")
    writeLines(c("@r1", "GGGAAAA", "+", "IIIIIII", "@r2", "GGGAAAA", "+", "IIIIIII"), fq)
    dset = count_fastq(fq, samples, constructs, read_structure="3B4M", verbose=FALSE)
    plot_read_count(dset)
}
```
