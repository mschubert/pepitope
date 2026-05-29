# Aggregate fusion VCFs into a table

Aggregate fusion VCFs into a table

## Usage

``` r
annotate_fusions(vr, txdb, asm)
```

## Arguments

- vr:

  A VRanges object with RNA fusions from readVcfAsRanges

- txdb:

  A transcription database, eg. AnnotationHub()\[\["AH100643"\]\]

- asm:

  A Genome sequence package object, eg. ::BSgenome.Hsapiens.NCBI.GRCh38

## Value

A DataFrame objects with fusions

## Examples

``` r
if (interactive()) {
    vcf = system.file("my_fusions.vcf", package="pepitope")
    vr = readVcfAsVRanges(vcf)
    txdb = AnnotationHub::AnnotationHub()[["AH100643"]]
    asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    annotate_fusions(vr, txdb, asm)
}
```
