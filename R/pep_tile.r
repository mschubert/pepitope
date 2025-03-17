#' Tile cDNA into peptide sequences
#'
#' @param df         A data.frame with variants to be tiled
#' @param tile_size  Oligo tiling size
#' @param tile_ov    Oligo tiling overlap
#'
#' @importFrom dplyr mutate filter select group_by ungroup rowwise n
#' @export
pep_tile = function(df, tile_size=93, tile_ov=45) {
    req = c("var_id", "mut_id", "gene_name", "gene_id", "tx_id", "cDNA")
    if (!all(req %in% colnames(df)))
        stop("Required column(s) not found: ", paste(setdiff(req, colnames(df)), collapse=", "))

    # tile peptides to have max `tile_size` nt length
    tile_cDNA = function(p) {
        ntile = ceiling((nchar(p)-tile_ov) / (tile_size-tile_ov))
        starts = round(seq(0, nchar(p)-tile_size, length.out=ntile)/3) * 3 + 1
        stopifnot((starts[-length(starts)]+tile_size - starts[-1]) >= tile_ov)
        lapply(starts, function(s) substr(p, s, s+tile_size-1))
    }

    pep = rowwise(df) %>%
            mutate(tiled = list(tile_cDNA(cDNA))) %>%
        ungroup() %>%
        mutate(n_tiles = sapply(tiled, length)) %>%
        tidyr::unnest(tiled) %>%
        mutate(tiled = unlist(tiled, use.names=FALSE),
               nt = nchar(tiled),
               peptide = as.character(translate(DNAStringSet(tiled), no.init.codon=TRUE))) %>%
        group_by(pep_id) %>%
            filter(!duplicated(tiled)) %>%
            mutate(pep_id = if(n()>1) paste(pep_id, seq_along(pep_id), sep="-") else pep_id) %>%
        ungroup()

#    stopifnot(pep$type[duplicated(pep$tiled)] == "ref")
    stopifnot(!any(duplicated(pep$pep_id)))
    stopifnot(vcountPattern("*", pep$peptide) == 0)
    pep
}
