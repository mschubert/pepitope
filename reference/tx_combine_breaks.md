# Combine break info from each possible left and right side transcript

Combine break info from each possible left and right side transcript

## Usage

``` r
tx_combine_breaks(vr, left, right)
```

## Arguments

- vr:

  A VRanges object with RNA fusions from \`readVcfAsVRanges\`

- left:

  List of DataFrame objects containing the 5' of the fusion

- right:

  List of DataFrame objects containing the 3' of the fusion

## Value

A DataFrame with fusion coordinates and sequence information
