pepitope: extract *pep*tide ep*itope*s
======================================

R package to extract the peptide context (flanking region around mutations)
from reference genome and variant (VCF) file.



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
* Save a report as _xlsx_ file

In code, this looks like the following:

```r
library(pepitope)

# genome and annotation
ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38

# example VCF file from the VariantAnnotation package
vr = readVcfAsVRanges("my.vcf.gz") |>
    filter_variants(min_cov=2, min_af=0.05, pass=TRUE)

# annotate and create report
ann = annotate_coding(vr, ens106, asm)
subs = subset_context(ann, ctx_codons=15)

report = make_report(ann, subs)
save_xlsx(report, "my_variants.xlsx")
```
