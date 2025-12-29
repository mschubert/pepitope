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
