# Minigene report

This vignette describes how to go from variants (and optionally gene
fusions) to the Minigene Library report.

We will usually work with standard processing pipelines for variant
calling, RNA-seq quantification, and RNA fusion calling. The best-known
collection of these workflows is [*NF-core*](https://nf-co.re/).

It provides the following pipelines for these tasks:

- [`sarek`](https://nf-co.re/sarek/) for variant calling, usually using
  *HaplotypeCaller* and *Mutect2*
- [`rnaseq`](https://nf-co.re/rnaseq/) for RNA-seq quantification,
  usually using *STAR* and *Salmon*
- [`rnafusion`](https://nf-co.re/rnafusion/) for RNA fusion
  quantification using multiple algorithms

Running these will usually require access to high performance computing
facilities, and you already get the processed files from your core
facility or sequencing provider.

Here, we are using example files provided with the package:

``` r
library(pepitope)

# SNPs and small indels, e.g. from 'sarek' nf-core pipeline
variant_vcf_file = system.file("my_variants.vcf", package="pepitope")

# Combined fusion VCF file, e.g. from 'rnafusion' nf-core pipeline
fusion_vcf_file = system.file("my_fusions.vcf", package="pepitope")
```

## Preparation

### Selecting the right reference genome

The variant calling and RNA-seq counts were mapped to a reference genome
and gene annotations. It is important to keep these consistent between
the NF-core processing pipelines and the Minigene Library annotation.

We are usually working with `BSgenome` (for the reference genome) and
`EnsDb` (for the gene annotations) objects. For human data, the most
widely used reference genome is `GRCh38`, and a recent *Ensembl*
annotation release.

With our test data, we know that `GRCh38` is the correct reference
genome and *Ensembl 106* is the correct version of gene annotations. We
can get both objects from the `BSgenome` and `AnnotationHub`
Bioconductor packages, respectively.

``` r
ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
#> loading from cache
#> require("ensembldb")
asm = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
```

A caveat here is that the chromosome prefixes need to be consistent
between the variants in the VCF file and the genome/gene annotations.
There are two “styles”, either UCSC (includes “chr” prefix) or
NCBI/Ensembl (without “chr” prefix).

The `sarek` pipline uses UCSC prefixes on the GRCh38 genome, so we need
to switch the gene annotation styles:

``` r
seqlevelsStyle(ens106) = "UCSC"
```

The correct styles will depend on how your VCF files were generated.

### Adding RNA expression

The [*NF-core* `rnaseq` workflow](https://nf-co.re/rnaseq/) will provide
two gene expression files, one for raw read counts and one for
transcripts per million (TPM). These contain all samples in a run, so we
need to subset them to the current sample we are interested in. These
files are usually called:

- `salmon.merged.gene_counts.tsv`
- `salmon.merged.gene_tpm.tsv`

We can combine and subset them the following way:

``` r
# note that this is not run in this example because we don't have RNA-seq data
counts = readr::read_tsv("salmon.merged.gene_counts.tsv") |>
    dplyr::select(gene_id, gene_name, count=SAMPLE)
tpm = readr::read_tsv("salmon.merged.gene_tpm.tsv") |>
    dplyr::select(gene_id, gene_name, tpm=SAMPLE)

rna_sample = inner_join(counts, tpm)
```

## SNPs and small indels

### Reading and filtering mutations

``` r
vr1 = readVcfAsVRanges(variant_vcf_file) |>
    filter_variants(min_cov=2, min_af=0.05, pass=TRUE)
```

Here, we can use the following filters:

- `min_cov = 2` – a variant needs to be covered by at least 2 reads
- `min_af = 0.05` – a variant needs to occur in at leat 5% of reads
- `pass = TRUE` – a variant needs to pass the standard QC filters
- `sample = SAMPLE` – only annotate the sample `SAMPLE`
- `chrs = "default"` – only use the common chromosomes for this organism

The resulting `vr1` object looks like the following:

### Annotating and subsetting expressed variants

``` r
ann = annotate_coding(vr1, ens106, asm)
subs = ann |>
#    filter_expressed(rna_sample, min_reads=1, min_tpm=0) |>
    subset_context(15)
```

Here, we are using the following filters:

- `min_reads = 1` – the gene needs to have at least one RNA read
- `min_tpm = 0` – we do not apply an additional TPM filter

In addition, we set the region of interest (context) to 15 codons up-
and downstream of the variant. Hence, a SNP will have a total length of
93 nucleotides (15\*3 + the SNP codon itself). An insertion will have
the inserted sequence and 15 codons, a deletion only 15 codons both
sides. A frameshift will have 15 codons upstream and the entire sequence
downstream until a STOP codon is reached. The latter may extend into the
3’ UTR.

The `subs` dataframe looks like the following:

## Fusion genes from RNA-seq

### Reading a fusion VCF

First we want to read the fusion genes from a combined `vcf` file like
the one produced by the [`rnafusion` NF-core
pipeline](https://nf-co.re/rnafusion/):

``` r
vr2 = readVcfAsVRanges(fusion_vcf_file) |>
    filter_fusions(min_reads=2, min_split_reads=1, min_tools=1)

seqlevelsStyle(vr2) = "UCSC"
```

Here, we are using the following filters:

- `min_reads = 2` – the fusion needs to be supported by 2 split or pair
  distance reads
- `min_split_reads = 1` – the fusion needs to be supported by at least
  one split read
- `min_tools = 1` – the fusion needs to be reported by at least one tool

### Annotating fusion genes

Next, we subset the peptide context analogous to the SNPs:

``` r
fus = annotate_fusions(vr2, ens106, asm) |>
    subset_context_fusion(15)
```

The `fus` table then looks like the following:

## Generating the Minigene Library

### Tiling cDNAs of interest into smaller peptides

``` r
tiled = make_peptides(subs, fus) |>
    pep_tile() |>
    remove_cutsite(BbsI="GAAGAC")
```

The `tiled` peptide table looks like the following:

### Saving the report file

We can then combine our generated tables into a report save it with
e.g. the `writexl` package. We can include all tables listed above:

``` r
report = make_report(ann, subs, fus, tiled)
writexl::write_xlsx(report, "report_file.xlsx")
```

This `.xlsx` file will contain the different tables as sheets. We will
use it as an annotation file in the quality control and screen steps.
