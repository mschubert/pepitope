# Make a variants report as named list of tables

Make a variants report as named list of tables

## Usage

``` r
make_report(vars, subs, fus = DataFrame(), tiled)
```

## Arguments

- vars:

  Variant results from \`annotate_coding()\`

- subs:

  Variants within mutation context from \`subset_context()\`

- fus:

  Variants within fusion context from \`subset_context_fusions()\`

- tiled:

  A data.frame of the tiled peptide sequences

## Value

A named list of data.frames for report output

## Examples

``` r
if (interactive()) {
    make_report(vars, subs, fus, tiled)
}
```
