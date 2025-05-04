pepitope: extract, qc and screen *pep*tide ep*itope*s
=====================================================

R package to extract the peptide context (flanking region around mutations)
from reference genome and variant (VCF) file.

The peptide context is defined as up- and downstream around the following:

* For SNPs, insertions, and deletions: all codons affected
* For frameshifts: the codon of the frameshift until the first stop codon

Installation
------------

The package is currently only available on Github, use the`remotes` package to
install:

```r
remotes::install_github("mschubert/pepitope")
```

Usage
-----

The _pepitope_ workflow usually consists of:

* Loading a genome and annotation, usually GRCh38 and Ensembl
* Loading a VCF variants file as `VRanges` object
* Annotate the protein-coding mutations
* Make a report of variants, coding changes, and tiled peptides

In code, this looks like the following:

```r
library(pepitope)

# genome and annotation
ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38

# read variants from VCF file and apply filters
vr = readVcfAsVRanges("my_variants.vcf.gz") |>
    filter_variants(min_cov=2, min_af=0.05, pass=TRUE)

# annotate and create report
ann = annotate_coding(vr, ens106, asm)

report = make_report(ann, ctx_codons=15)
# writexl::write_xlsx(report, "my_variants.xlsx")
```
