#' Tile cDNA into peptide sequences
#'
#' @param subs       A GRanges object of context-subset protein-coding variants
#' @param tile_size  Oligo tiling size
#' @param tile_ov    Oligo tiling overlap
#'
#' @importFrom dplyr mutate filter select group_by ungroup rowwise n
#' @importFrom S4Vectors mcols
#' @export
pep_tile = function(subs, tile_size=93, tile_ov=45) {
    req = c("var_id", "mut_id", "gene_name", "GENEID", "tx_name", "ref_nuc", "alt_nuc")
    if (!all(req %in% colnames(mcols(subs))))
        stop("Required column(s) not found: ", paste(setdiff(req, colnames(df)), collapse=", "))

    df = as.data.frame(subs) %>%
        select(var_id, mut_id, gene_name, gene_id=GENEID, tx_id=tx_name,
               ref=ref_nuc, alt=alt_nuc) %>%
        tidyr::pivot_longer(c(ref, alt), names_to="type", values_to="cDNA") %>%
        mutate(pep_id = ifelse(type == "alt", mut_id, sub("([0-9]+)[a-zA-Z*]+$", "\\1", mut_id)))

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
