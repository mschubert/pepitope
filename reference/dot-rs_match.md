# Match inferred barcode segments in individual reads

Match inferred barcode segments in individual reads

## Usage

``` r
.rs_match(reads, structure, sample_barcodes, construct_barcodes)
```

## Arguments

- reads:

  Read sequences as a character vector or \`DNAStringSet\`

- structure:

  An inferred read structure string, including \`"unknown"\`

- sample_barcodes:

  Known sample barcodes

- construct_barcodes:

  Known construct barcodes

## Value

A \`data.frame\` with read indices, barcode types, ranges, orientations,
and logical \`present\` values; unknown structures return zero rows
