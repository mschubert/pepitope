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
