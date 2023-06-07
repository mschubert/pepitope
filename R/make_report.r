#' Make a variants report as named list of tables
#'
#' @param res   Variant results
#' @param subs  Filtered variant results
#' @param fus   Fusion results
#' @param ...   Parameters passed to `pep_tile` (eg. tile_size, tile_ov)
#' @export
make_report = function(res, subs, fus, ...) {
    gr2df = function(gr) as_tibble(as.data.frame(gr)) %>%
        select(var_id, everything()) %>%
        arrange(order(gtools::mixedorder(var_id)))

    pep1 = pep_tile(subs, ...)
#    pep2 = pep_tile(fus, ...)

    list(`All Variants` = res %>% select(-CDSLOC) %>% gr2df() %>%
            dplyr::select(-tx_name, -(ref_nuc:alt_prot)) %>% distinct(),
         `Unique Protein-Coding` = gr2df(subs) %>% select(var_id, mut_id, everything()),
         `93 nt Peptides` = pep %>% select(var_id, mut_id, pep_id,
            gene_id:cDNA, n_tiles, BbsI_replaced, tiled, nt, peptide))
}
