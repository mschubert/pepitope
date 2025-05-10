pepitope: extract, qc and screen *pep*tide ep*itope*s
=====================================================

This R package is used to:

* extract peptides with flanking region around mutations ([vignette](https://mschubert.github.io/pepitope/articles/variant.html))
* run quality control on sequencing of these libraries ([vignette](https://mschubert.github.io/pepitope/articles/qc.html))
* perform differential abundance testing of co-culture screens ([vignette](https://mschubert.github.io/pepitope/articles/screen.html))

Installation
------------

We require `R>=4.5.0`, as this includes important fixes on Bioconductor.

The package is currently only available on Github, use the`remotes` package to
install:

```r
remotes::install_github("mschubert/pepitope")
```

In addition, we need the [Rust](https://www.rust-lang.org/tools/install)
programming language and its `cargo` package manager, which we use to install
[`fqtk`](https://github.com/fulcrumgenomics/fqtk) and
[`guide-counter`](https://github.com/fulcrumgenomics/guide-counter):

```sh
cargo install fqtk guide-counter
```

We can check if the installation works by running, in R:

```r
library(pepitope)
Sys.which(c("fqtk", "guide-counter")) # should print paths, not ""
```

Usage
-----

#### Generating peptide constructs ([vignette](https://mschubert.github.io/pepitope/articles/variant.html))

* Load a genome and annotation, usually GRCh38 and Ensembl
* Load a VCF variants file as `VRanges` object
* Annotate the protein-coding mutations
* Optionally, loading a fusion VCF and annotating those
* Subsetting the peptide context around each mutation
* Make a report of variants, coding changes, and tiled peptides

<details><summary><b>Code example</b></summary>

```r
library(pepitope)

# genome and annotation
ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(ens106) = "UCSC"

# read variants from VCF file, apply filters and annotate
variant_vcf_file = system.file("my_variants.vcf", package="pepitope")
vr = readVcfAsVRanges(variant_vcf_file) |>
    filter_variants(min_cov=2, min_af=0.05, pass=TRUE)
ann = annotate_coding(vr, ens106, asm)
subs = ann |>
#    filter_expressed(rna_sample, min_reads=1, min_tpm=0) |>
    subset_context(15)

# read fusion variants, apply filters and annotate
fusion_vcf_file = system.file("my_fusions.vcf", package="pepitope")
vr2 = readVcfAsVRanges(fusion_vcf_file) |>
    filter_fusions(min_reads=2, min_split_reads=1, min_tools=1)
seqlevelsStyle(vr2) = "UCSC"
fus = annotate_fusions(vr2, ens106, asm) |>
    subset_context_fusion(15)

# create construct tables and make a report
tiled = make_peptides(subs, fus) |>
    pep_tile() |>
    remove_cutsite(BbsI="GAAGAC")

report = make_report(ann, subs, fus, tiled)
writexl::write_xlsx(report, "my_variants.xlsx")
```

</details>

#### Creating construct library (wetlab)

* Decide on barcodes for each of those constructs
* Add them in the annotation sheets as `barcode` or `barcode_{1,2}` etc.
* Order these constructs as gene blocks and clone them into expression vectors
* Transduce target cells with this peptide construct library

<details><summary><b>Code example</b></summary>

Normally, we want to create an `xlsx` document with all libraries in use:

```r
# this file is manually created from the output of step 1
fname = "my_combined_barcoded_file.xlsx"
sheets = readxl::excel_sheets(fname)
all_constructs = sapply(sheets, readxl::read_xlsx, path=fname, simplify=FALSE)
```

Here, we use our example library instead:

```r
# creating barcoded constructs
lib = "https://raw.githubusercontent.com/hawkjo/freebarcodes/master/barcodes/barcodes12-1.txt"
valid_barcodes = readr::read_tsv(lib, col_names=FALSE)$X1
all_constructs = example_peptides(valid_barcodes)
```

</details>

#### Performing quality control on construct library sequencing ([vignette](https://mschubert.github.io/pepitope/articles/qc.html))

* Check the quality of the construct libraries
* Check the quality of the transduced target cells
* Check the quality of the co-culture screens

<details><summary><b>Code example</b></summary>

Quality Control can be performed on the library, the target cells, or
the co-culture screen:

```r
# demultiplexing and counting example data
sample_sheet = system.file("my_samples.tsv", package="pepitope")
fastq_file = example_fastq(sample_sheet, all_constructs)
temp_dir = demux_fq(fastq_file, sample_sheet, read_structures="7B+T")
dset = count_bc(temp_dir, all_constructs, valid_barcodes)

# quality control plots
plot_barcode_overlap(all_constructs, valid_barcodes)
plot_reads(dset)
plot_distr(dset)
```

</details>

#### Differential abundance testing of co-culture screens ([vignette](https://mschubert.github.io/pepitope/articles/screen.html))

<details><summary><b>Code example</b></summary>

Differential abundance testing works the following way:

```r
# perform abundance testing and plot results
res = screen_calc(dset, list(c("Sample", "Mock")))
plot_screen(res$`Sample vs Mock`)
```

</details>
