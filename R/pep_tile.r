#' Tile cDNA into peptide sequences
#'
#' @param subs       A GRanges object
#' @param tile_size  Oligo tiling size
#' @param tile_ov    Oligo tiling overlap
#' @export
pep_tile = function(subs, tile_size=93, tile_ov=45) {
    req = c("var_id", "mut_id", "gene_name", "GENEID", "tx_name", "ref_nuc", "alt_nuc")
    if (!all(req %in% colnames(mcols(subs))))
        stop("Required column(s) not found: ", paste(setdiff(req, colnames(subs)), collapse=", "))

    # tile peptides to have max `tile_size` nt length
    tile_cDNA = function(p) {
        ntile = ceiling((nchar(p)-tile_ov) / (tile_size-tile_ov))
        starts = round(seq(0, nchar(p)-tile_size, length.out=ntile)/3) * 3 + 1
        stopifnot((starts[-length(starts)]+tile_size - starts[-1]) >= tile_ov)
        lapply(starts, function(s) substr(p, s, s+tile_size-1))
    }
    pep = with(subs, tibble(var_id=subs$var_id, mut_id=subs$mut_id,
                            gene_name=subs$gene_name, gene_id=subs$GENEID,
                            tx_id=subs$tx_name, ref=as.character(subs$ref_nuc),
                            alt=as.character(subs$alt_nuc))) %>%
        tidyr::pivot_longer(c(ref, alt), names_to="type", values_to="cDNA") %>%
        rowwise() %>% mutate(tiled = list(tile_cDNA(cDNA))) %>% ungroup() %>%
        mutate(n_tiles = sapply(tiled, length)) %>%
        tidyr::unnest(tiled) %>%
        mutate(pep_id = ifelse(type == "alt", mut_id, sub("([0-9]+)[a-zA-Z]+$", "\\1", mut_id)),
               tiled = unlist(tiled, use.names=FALSE),
               nt = nchar(tiled),
               peptide = as.character(translate(DNAStringSet(tiled), no.init.codon=TRUE))) %>%
        group_by(pep_id) %>%
            mutate(pep_id = if(n()>1) paste(pep_id, seq_along(pep_id), sep="-") else pep_id) %>%
        ungroup()

#    stopifnot(pep$type[duplicated(pep$tiled)] == "ref")
    stopifnot(!any(duplicated(pep$pep_id)))
    pep = pep[!duplicated(pep$tiled),]

    pep$BbsI_replaced = vcountPattern("GAAGAC", pep$tiled) + vcountPattern("GTCTTC", pep$tiled)
    pep$tiled = sapply(pep$tiled, remove_cutsite, site="GAAGAC", seed=178529, USE.NAMES=FALSE)
#    stopifnot(pep$peptide == as.character(translate(DNAStringSet(pep$tiled), no.init.codon=TRUE)))

    pep
}
