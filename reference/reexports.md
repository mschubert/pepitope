# Objects exported from other packages

These objects are imported from other packages. Follow the links below
to see their documentation.

## Usage

``` r
assay(x, i, withDimnames = TRUE, ...)

colData(x, ...)

rowData(x, use.names = TRUE, ...)

genome(x)

seqinfo(x)

seqlevels(x)

seqlevelsStyle(x)

seqlevelsStyle(x) <- value

readVcfAsVRanges(x, genome, param = ScanVcfParam(), use.names = FALSE, ...)
```

## Details

- GenomeInfoDb:

  [`genome`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html),
  [`seqinfo`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html),
  [`seqlevels`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html),
  [`seqlevelsStyle()`](https://rdrr.io/pkg/GenomeInfoDb/man/seqlevelsStyle.html),
  [`seqlevelsStyle<-()`](https://rdrr.io/pkg/GenomeInfoDb/man/seqlevelsStyle.html)

- SummarizedExperiment:

  [`assay()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
  [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
  [`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

- VariantAnnotation:

  [`readVcfAsVRanges()`](https://rdrr.io/pkg/VariantAnnotation/man/VRanges-class.html)
