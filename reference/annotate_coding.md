# Annotate VCF variants with coding changes

Annotate VCF variants with coding changes

## Usage

``` r
annotate_coding(vr, txdb, asm)
```

## Arguments

- vr:

  A VRanges object with SNVs and small indels

- txdb:

  TxDb or EnsDb object

- asm:

  Genomic sequence BSGenome object

## Value

A GRanges object with annotated variants

## Examples

``` r
if (interactive()) {
    vcf = system.file("my_variants.vcf", package="pepitope")
    vr = readVcfAsVRanges(vcf)
    txdb = AnnotationHub::AnnotationHub()[["AH100643"]]
    asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    annotate_coding(vr, txdb, asm)
}
```
