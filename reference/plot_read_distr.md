# Plot the read distribution across barcodes

Plot the read distribution across barcodes

## Usage

``` r
plot_read_distr(dset)

plot_distr(dset)
```

## Arguments

- dset:

  The \`SummarizedExperiment\` object from \`count_fastq\`

## Value

A \`ggplot2\` object with cumulative read distribution plots

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
    plot_read_distr(dset)
}
```
