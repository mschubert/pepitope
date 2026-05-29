# Filter a fusion VRanges object by number of reads and tools

Filter a fusion VRanges object by number of reads and tools

## Usage

``` r
filter_fusions(vr, min_reads = NULL, min_split_reads = NULL, min_tools = NULL)
```

## Arguments

- vr:

  A VRanges object with RNA fusions from readVcfAsRanges

- min_reads:

  The minimum number of linked read support for a fusion

- min_split_reads:

  The minimum number of split read support for a fusion

- min_tools:

  The minimum number of tools that identify a fusion

## Value

A filtered VRanges object

## Examples

``` r
vr = VariantAnnotation::VRanges("chr1", IRanges::IRanges(1, 1), ref="A", alt="T")
vr$DV = 3L
vr$RV = 2L
vr$TOOL_HITS = IRanges::IntegerList(2L)
filter_fusions(vr, min_reads=4, min_split_reads=2, min_tools=1)
#> VRanges object with 1 range and 3 metadata columns:
#>       seqnames    ranges strand         ref              alt     totalDepth
#>          <Rle> <IRanges>  <Rle> <character> <characterOrRle> <integerOrRle>
#>   [1]     chr1         1      *           A                T           <NA>
#>             refDepth       altDepth   sampleNames softFilterMatrix |        DV
#>       <integerOrRle> <integerOrRle> <factorOrRle>         <matrix> | <integer>
#>   [1]           <NA>           <NA>          <NA>                  |         3
#>              RV     TOOL_HITS
#>       <integer> <IntegerList>
#>   [1]         2             2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#>   hardFilters: NULL
```
