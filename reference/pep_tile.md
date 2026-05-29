# Tile cDNA into peptide sequences

Tile cDNA into peptide sequences

## Usage

``` r
pep_tile(peptides, tile_size = 93, tile_ov = 45)
```

## Arguments

- peptides:

  A \`data.frame\` with context-subset peptide/minigene data

- tile_size:

  Oligo tiling size

- tile_ov:

  Oligo tiling overlap

## Value

A data.frame with tiled nucleotide and peptide sequences

## Examples

``` r
peptides = data.frame(var_id="v1", mut_id="M1_A1V", gene_name="M1",
    gene_id="gene1", tx_id="tx1", pep_id="M1_A1V", cDNA="ATGGCCGCCGCC")
pep_tile(peptides, tile_size=9, tile_ov=3)
#> # A tibble: 2 × 10
#>   var_id mut_id gene_name gene_id tx_id pep_id   tiled     n_tiles    nt peptide
#>   <chr>  <chr>  <chr>     <chr>   <chr> <chr>    <chr>       <int> <int> <chr>  
#> 1 v1     M1_A1V M1        gene1   tx1   M1_A1V-1 ATGGCCGCC       2     9 MAA    
#> 2 v1     M1_A1V M1        gene1   tx1   M1_A1V-2 GCCGCCGCC       2     9 AAA    
```
