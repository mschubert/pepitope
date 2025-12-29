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
