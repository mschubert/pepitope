# Get a list of transcripts with their CDS exons overlapping the break

Get a list of transcripts with their CDS exons overlapping the break

## Usage

``` r
cds_by_break(gr, txdb, cds, type = "left")
```

## Arguments

- gr:

  GenomicRanges object of break location

- txdb:

  A transcription database, eg. AnnotationHub()\[\["AH100643"\]\]

- cds:

  A list of exon coordinates by gene from \`cdsBy(txdb)\`

- type:

  Whether we want info for the 'left' or 'right' side of the break

## Value

A named list of transcript GRanges objects with CDS exons
