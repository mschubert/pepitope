# Adds sequence information to break transcripts

Adds sequence information to break transcripts

## Usage

``` r
add_seq_info(gr, cds_break, asm, txdb, tx)
```

## Arguments

- gr:

  GenomicRanges object of break location

- cds_break:

  A list of transcripts overlapping break from \`cds_by_break\`

- asm:

  A Genome sequence package object, eg. ::BSgenome.Hsapiens.NCBI.GRCh38

- txdb:

  A transcription database, eg. AnnotationHub()\[\["AH100643"\]\]

- tx:

  A list of transcripts obtained from \`transcripts(txdb)\`

## Value

A DataFrame with sequence information
