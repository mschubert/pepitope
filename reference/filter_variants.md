# Make results report to save as xlsx sheets (full, filtered, peptides)

Make results report to save as xlsx sheets (full, filtered, peptides)

## Usage

``` r
filter_variants(
  vr,
  ...,
  min_cov = 2,
  min_af = 0.05,
  pass = TRUE,
  sample = NULL,
  chrs = NULL
)
```

## Arguments

- vr:

  A VRanges object from \`readVcfAsVRanges\`

- ...:

  Force filters by name (ignored)

- min_cov:

  Minimum number of reads to span the ALT allele

- min_af:

  Minimum allele frequency of the ALT allele

- pass:

  Whether to only include softFilterMatrix PASS

- sample:

  Only include if in \`sampleNames(vr)\` (required if more than one
  present)

- chrs:

  Either "default" or a character vector of chromosome names

## Value

A filtered VRanges object

## Examples

``` r
vcf = system.file("my_variants.vcf", package="pepitope")
vr = readVcfAsVRanges(vcf)
filter_variants(vr, min_cov=2, min_af=0.05, pass=TRUE)
#> VRanges object with 32 ranges and 1 metadata column:
#>        seqnames    ranges strand         ref              alt     totalDepth
#>           <Rle> <IRanges>  <Rle> <character> <characterOrRle> <integerOrRle>
#>    [1]     chr1 114713908      *           T                A            100
#>    [2]     chr1 114713908      *           T                C            100
#>    [3]     chr2 177234082      *           G                A            100
#>    [4]     chr2 198267522      *           A                G            100
#>    [5]     chr2 208248388      *           G                A            100
#>    ...      ...       ...    ...         ...              ...            ...
#>   [28]    chr15  90088702      *           C                T            100
#>   [29]    chr17   3123250      *           A                T            100
#>   [30]    chr17   7675087      *           C                T            100
#>   [31]    chr17  43081873      *           C               CC            100
#>   [32]    chr18  51065548      *           C                T            100
#>              refDepth       altDepth   sampleNames softFilterMatrix |      QUAL
#>        <integerOrRle> <integerOrRle> <factorOrRle>         <matrix> | <numeric>
#>    [1]             90             10        SAMPLE                  |        NA
#>    [2]             90             10        SAMPLE                  |        NA
#>    [3]             90             10        SAMPLE                  |        NA
#>    [4]             90             10        SAMPLE                  |        NA
#>    [5]             90             10        SAMPLE                  |        NA
#>    ...            ...            ...           ...              ... .       ...
#>   [28]             90             10        SAMPLE                  |        NA
#>   [29]             90             10        SAMPLE                  |        NA
#>   [30]             90             10        SAMPLE                  |        NA
#>   [31]             90             10        SAMPLE                  |        NA
#>   [32]             90             10        SAMPLE                  |        NA
#>   -------
#>   seqinfo: 14 sequences from an unspecified genome; no seqlengths
#>   hardFilters: NULL
```
