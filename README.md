pepitope: extract, qc and screen *pep*tide ep*itope*s
=====================================================

This R package is used to:

* [**extract peptides**](#generating-peptide-constructs) with flanking region around mutations
* run [**quality control**](#performing-quality-control-on-construct-library-sequencing) on sequencing of these libraries
* perform [**differential abundance**](#differential-abundance-testing-of-co-culture-screens) testing of co-culture screens

Installation
------------

We require `R>=4.5.0`, as this includes important fixes on Bioconductor.

The package is currently only available on Github, use the`remotes` package to
install:

```r
# run this in R
if (!require(remotes))
    install.packages("remotes")
remotes::install_github("mschubert/pepitope")
```

In addition, we need the [Rust](https://www.rust-lang.org/tools/install)
programming language and its `cargo` package manager, which we use to install
[`fqtk`](https://github.com/fulcrumgenomics/fqtk) and
[`guide-counter`](https://github.com/fulcrumgenomics/guide-counter):

```sh
# run this in your terminal
cargo install fqtk guide-counter
```

We can check if the installation works by running:

```r
# run this in R
library(pepitope)
Sys.which(c("fqtk", "guide-counter")) # should print paths, not ""
```

Usage
-----

### Generating peptide constructs

Here we have sequenced the DNA (and optionally RNA) of a patient and
identified the variants in a `.vcf` file. We now want to extract the
reference and mutated alternative sequences including their flanking
regions into a summary report. The steps are:

* Load a genome and annotation, usually GRCh38 and Ensembl
* Load a VCF variants file as `VRanges` object and annotate the protein-coding mutations
* Optionally, load a fusion VCF and annotating those
* Subset the peptide context around each mutation
* Make a report of variants, coding changes, and tiled peptides

More information can be found in the [*Variant calling* vignette 🔗](https://mschubert.github.io/pepitope/articles/variant.html).

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

### Creating construct library (wetlab)

We want to express the sequences (minigenes) including their flanking regions
(context) in target cells that will be used in a co-culture screen with T-cells.
For this, we first need to add a barcode to each construct and then order
them as gene blocks and transduce them into the target cells. The steps are:

* Decide on barcodes for each of those constructs
* Add them in the annotation sheets as `barcode` or `barcode_{1,2}` etc.
* Order these constructs as gene blocks and clone them into expression vectors
* Transduce target cells with this peptide construct library

<details><summary><b>Code example</b></summary>

Normally, we want to create an `xlsx` document where each sheet is one librar we use:

```r
# this file is manually created from the output of step 1
fname = "my_combined_barcoded_file.xlsx"
sheets = readxl::excel_sheets(fname)
all_constructs = sapply(sheets, readxl::read_xlsx, path=fname, simplify=FALSE)
```

Here, we use our example libraries instead:

```r
# creating barcoded constructs
lib = "https://raw.githubusercontent.com/hawkjo/freebarcodes/master/barcodes/barcodes12-1.txt"
valid_barcodes = readr::read_tsv(lib, col_names=FALSE)$X1
all_constructs = example_peptides(valid_barcodes)
```

</details>

### Performing quality control on construct library sequencing

In each step of generating the target cells expressing the reference and mutated
versions of each peptide, we want to make sure our library is well-represented.
For this, we will check if all constructs that should be in there are, and whether
they are present in a similar enough amount. The steps are:

* Check the quality of the construct libraries
* Check the quality of the transduced target cells
* Check the quality of the co-culture screens

More information can be found in the [*QC* vignette 🔗](https://mschubert.github.io/pepitope/articles/qc.html)

<details><summary><b>Code example</b></summary>

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

### Differential abundance testing of co-culture screens

Finally, we co-culture our target cells with T-cells expressing a variety of
TCRs with our expressed peptide libraries to find the reactive ones. Those will
be visible by decreasing in abundance more than the reference peptides compared
to a mock-transduced population that was cultured the same way. The steps are:

* Calculate the differential abundance of peptide barcodes
* Plot the results to identify peptides recognized by T-cells

More information can be found in the [*screen* vignette 🔗](https://mschubert.github.io/pepitope/articles/screen.html).

<details><summary><b>Code example</b></summary>

```r
# perform abundance testing and plot results
res = screen_calc(dset, list(c("Sample", "Mock")))
plot_screen(res$`Sample vs Mock`)
```

</details>
