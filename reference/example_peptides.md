# Create example peptide sheets for multiple samples

Create example peptide sheets for multiple samples

## Usage

``` r
example_peptides(valid_barcodes)
```

## Arguments

- valid_barcodes:

  A character vector of valid barcodes

## Value

A named list of peptide/minigene constructs with barcodes

## Examples

``` r
valid_barcodes = example_barcodes(1000)
example_peptides(valid_barcodes)
#> $pat1
#> # A tibble: 67 × 13
#>    var_id             mut_id       pep_id        pep_type gene_name
#>    <chr>              <chr>        <chr>         <chr>    <chr>    
#>  1 chr1:114713908_T/A NRAS_Q61     NRAS_Q61      ref      NRAS     
#>  2 chr1:114713908_T/A NRAS_Q61     NRAS_Q61L     alt      NRAS     
#>  3 chr1:114713908_T/C NRAS_Q61     NRAS_Q61R     alt      NRAS     
#>  4 chr2:208248388_G/A IDH1_R132    IDH1_R132     ref      IDH1     
#>  5 chr2:208248388_G/A IDH1_R132    IDH1_R132L    alt      IDH1     
#>  6 chr3:179234297_A/G PIK3CA_H1047 PIK3CA_H1047  ref      PIK3CA   
#>  7 chr3:179234297_A/G PIK3CA_H1047 PIK3CA_H1047R alt      PIK3CA   
#>  8 chr4:54285926_A/T  PDGFRA_D842  PDGFRA_D842   ref      PDGFRA   
#>  9 chr4:54285926_A/T  PDGFRA_D842  PDGFRA_D842V  alt      PDGFRA   
#> 10 chr4:54733155_A/T  KIT_D816     KIT_D816      ref      KIT      
#>    gene_id         tx_id         tiled n_tiles    nt peptide barcode_1 barcode_2
#>    <chr>           <chr>         <chr>   <dbl> <dbl> <chr>   <chr>     <chr>    
#>  1 ENSG00000213281 ENST00000369… ATAG…       1    93 IDGETC… AAAAAAAA… TAACAAAA…
#>  2 ENSG00000213281 ENST00000369… ATAG…       1    93 IDGETC… CAAAAAAA… ACACAAAA…
#>  3 ENSG00000213281 ENST00000369… ATAG…       1    93 IDGETC… GAAAAAAA… CCACAAAA…
#>  4 ENSG00000138413 ENST00000345… ATCC…       1    93 IPRLVS… TAAAAAAA… GCACAAAA…
#>  5 ENSG00000138413 ENST00000345… ATCC…       1    93 IPRLVS… ACAAAAAA… TCACAAAA…
#>  6 ENSG00000121879 ENST00000263… GAGC…       1    93 EQEALE… CCAAAAAA… AGACAAAA…
#>  7 ENSG00000121879 ENST00000263… GAGC…       1    93 EQEALE… GCAAAAAA… CGACAAAA…
#>  8 ENSG00000134853 ENST00000257… GCAC…       1    93 AQGKIV… TCAAAAAA… GGACAAAA…
#>  9 ENSG00000134853 ENST00000257… GCAC…       1    93 AQGKIV… AGAAAAAA… TGACAAAA…
#> 10 ENSG00000157404 ENST00000288… ACTC…       1    93 THGRIT… CGAAAAAA… ATACAAAA…
#> # ℹ 57 more rows
#> 
#> $pat2
#> # A tibble: 46 × 13
#>    var_id             mut_id       pep_id        pep_type gene_name
#>    <chr>              <chr>        <chr>         <chr>    <chr>    
#>  1 chr1:114713908_T/A NRAS_Q61     NRAS_Q61      ref      NRAS     
#>  2 chr1:114713908_T/A NRAS_Q61     NRAS_Q61L     alt      NRAS     
#>  3 chr1:114713908_T/C NRAS_Q61     NRAS_Q61R     alt      NRAS     
#>  4 chr3:179234297_A/G PIK3CA_H1047 PIK3CA_H1047  ref      PIK3CA   
#>  5 chr3:179234297_A/G PIK3CA_H1047 PIK3CA_H1047R alt      PIK3CA   
#>  6 chr4:54285926_A/T  PDGFRA_D842  PDGFRA_D842   ref      PDGFRA   
#>  7 chr4:54285926_A/T  PDGFRA_D842  PDGFRA_D842V  alt      PDGFRA   
#>  8 chr4:54733155_A/T  KIT_D816     KIT_D816      ref      KIT      
#>  9 chr4:54733155_A/T  KIT_D816     KIT_D816V     alt      KIT      
#> 10 chr4:152326137_C/T FBXW7_R505   FBXW7_R505    ref      FBXW7    
#>    gene_id         tx_id         tiled n_tiles    nt peptide barcode_1 barcode_2
#>    <chr>           <chr>         <chr>   <dbl> <dbl> <chr>   <chr>     <chr>    
#>  1 ENSG00000213281 ENST00000369… ATAG…       1    93 IDGETC… GCAGAAAA… ACTGAAAA…
#>  2 ENSG00000213281 ENST00000369… ATAG…       1    93 IDGETC… TCAGAAAA… CCTGAAAA…
#>  3 ENSG00000213281 ENST00000369… ATAG…       1    93 IDGETC… AGAGAAAA… GCTGAAAA…
#>  4 ENSG00000121879 ENST00000263… GAGC…       1    93 EQEALE… CGAGAAAA… TCTGAAAA…
#>  5 ENSG00000121879 ENST00000263… GAGC…       1    93 EQEALE… GGAGAAAA… AGTGAAAA…
#>  6 ENSG00000134853 ENST00000257… GCAC…       1    93 AQGKIV… TGAGAAAA… CGTGAAAA…
#>  7 ENSG00000134853 ENST00000257… GCAC…       1    93 AQGKIV… ATAGAAAA… GGTGAAAA…
#>  8 ENSG00000157404 ENST00000288… ACTC…       1    93 THGRIT… CTAGAAAA… TGTGAAAA…
#>  9 ENSG00000157404 ENST00000288… ACTC…       1    93 THGRIT… GTAGAAAA… ATTGAAAA…
#> 10 ENSG00000109670 ENST00000281… ACAG…       1    93 TGQCLH… TTAGAAAA… CTTGAAAA…
#> # ℹ 36 more rows
#> 
#> $pat3
#> # A tibble: 47 × 13
#>    var_id             mut_id       pep_id        pep_type gene_name
#>    <chr>              <chr>        <chr>         <chr>    <chr>    
#>  1 chr2:208248388_G/A IDH1_R132    IDH1_R132     ref      IDH1     
#>  2 chr2:208248388_G/A IDH1_R132    IDH1_R132L    alt      IDH1     
#>  3 chr3:179234297_A/G PIK3CA_H1047 PIK3CA_H1047  ref      PIK3CA   
#>  4 chr3:179234297_A/G PIK3CA_H1047 PIK3CA_H1047R alt      PIK3CA   
#>  5 chr4:54285926_A/T  PDGFRA_D842  PDGFRA_D842   ref      PDGFRA   
#>  6 chr4:54285926_A/T  PDGFRA_D842  PDGFRA_D842V  alt      PDGFRA   
#>  7 chr4:54733155_A/T  KIT_D816     KIT_D816      ref      KIT      
#>  8 chr4:54733155_A/T  KIT_D816     KIT_D816V     alt      KIT      
#>  9 chr4:152326137_C/T FBXW7_R505   FBXW7_R505    ref      FBXW7    
#> 10 chr4:152326137_C/T FBXW7_R505   FBXW7_R505S   alt      FBXW7    
#>    gene_id         tx_id         tiled n_tiles    nt peptide barcode_1 barcode_2
#>    <chr>           <chr>         <chr>   <dbl> <dbl> <chr>   <chr>     <chr>    
#>  1 ENSG00000138413 ENST00000345… ATCC…       1    93 IPRLVS… GAGTAAAA… CACACAAA…
#>  2 ENSG00000138413 ENST00000345… ATCC…       1    93 IPRLVS… TAGTAAAA… GACACAAA…
#>  3 ENSG00000121879 ENST00000263… GAGC…       1    93 EQEALE… ACGTAAAA… TACACAAA…
#>  4 ENSG00000121879 ENST00000263… GAGC…       1    93 EQEALE… CCGTAAAA… ACCACAAA…
#>  5 ENSG00000134853 ENST00000257… GCAC…       1    93 AQGKIV… GCGTAAAA… CCCACAAA…
#>  6 ENSG00000134853 ENST00000257… GCAC…       1    93 AQGKIV… TCGTAAAA… GCCACAAA…
#>  7 ENSG00000157404 ENST00000288… ACTC…       1    93 THGRIT… AGGTAAAA… TCCACAAA…
#>  8 ENSG00000157404 ENST00000288… ACTC…       1    93 THGRIT… CGGTAAAA… AGCACAAA…
#>  9 ENSG00000109670 ENST00000281… ACAG…       1    93 TGQCLH… GGGTAAAA… CGCACAAA…
#> 10 ENSG00000109670 ENST00000281… ACAG…       1    93 TGQCLH… TGGTAAAA… GGCACAAA…
#> # ℹ 37 more rows
#> 
#> $common
#> # A tibble: 46 × 13
#>    var_id             mut_id       pep_id        pep_type gene_name
#>    <chr>              <chr>        <chr>         <chr>    <chr>    
#>  1 chr1:114713908_T/A NRAS_Q61     NRAS_Q61      ref      NRAS     
#>  2 chr1:114713908_T/A NRAS_Q61     NRAS_Q61L     alt      NRAS     
#>  3 chr1:114713908_T/C NRAS_Q61     NRAS_Q61R     alt      NRAS     
#>  4 chr2:208248388_G/A IDH1_R132    IDH1_R132     ref      IDH1     
#>  5 chr2:208248388_G/A IDH1_R132    IDH1_R132L    alt      IDH1     
#>  6 chr3:179234297_A/G PIK3CA_H1047 PIK3CA_H1047  ref      PIK3CA   
#>  7 chr3:179234297_A/G PIK3CA_H1047 PIK3CA_H1047R alt      PIK3CA   
#>  8 chr4:54285926_A/T  PDGFRA_D842  PDGFRA_D842   ref      PDGFRA   
#>  9 chr4:54285926_A/T  PDGFRA_D842  PDGFRA_D842V  alt      PDGFRA   
#> 10 chr4:54733155_A/T  KIT_D816     KIT_D816      ref      KIT      
#>    gene_id         tx_id         tiled n_tiles    nt peptide barcode_1 barcode_2
#>    <chr>           <chr>         <chr>   <dbl> <dbl> <chr>   <chr>     <chr>    
#>  1 ENSG00000213281 ENST00000369… ATAG…       1    93 IDGETC… AAACCAAA… GTGCCAAA…
#>  2 ENSG00000213281 ENST00000369… ATAG…       1    93 IDGETC… CAACCAAA… TTGCCAAA…
#>  3 ENSG00000213281 ENST00000369… ATAG…       1    93 IDGETC… GAACCAAA… AATCCAAA…
#>  4 ENSG00000138413 ENST00000345… ATCC…       1    93 IPRLVS… TAACCAAA… CATCCAAA…
#>  5 ENSG00000138413 ENST00000345… ATCC…       1    93 IPRLVS… ACACCAAA… GATCCAAA…
#>  6 ENSG00000121879 ENST00000263… GAGC…       1    93 EQEALE… CCACCAAA… TATCCAAA…
#>  7 ENSG00000121879 ENST00000263… GAGC…       1    93 EQEALE… GCACCAAA… ACTCCAAA…
#>  8 ENSG00000134853 ENST00000257… GCAC…       1    93 AQGKIV… TCACCAAA… CCTCCAAA…
#>  9 ENSG00000134853 ENST00000257… GCAC…       1    93 AQGKIV… AGACCAAA… GCTCCAAA…
#> 10 ENSG00000157404 ENST00000288… ACTC…       1    93 THGRIT… CGACCAAA… TCTCCAAA…
#> # ℹ 36 more rows
#> 
```
