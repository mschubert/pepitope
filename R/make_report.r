#' Make a variants report as named list of tables
#'
#' @param vars   Variant results
#' @param subs  Filtered variant results
#' @param fus   Fusion results
#' @param ...   Parameters passed to `pep_tile` (eg. tile_size, tile_ov)
#'
#' @importFrom dplyr bind_rows
#' @export
make_report = function(vars, subs, fus=DataFrame(), ...) {
    gr2df = function(gr) as_tibble(as.data.frame(gr)) %>%
        select(var_id, everything()) %>%
        arrange(order(gtools::mixedorder(var_id)))

    pep = as.data.frame(subs) %>%
        select(var_id, mut_id, gene_name, gene_id=GENEID, tx_id=tx_name,
               ref=ref_nuc, alt=alt_nuc) %>%
        tidyr::pivot_longer(c(ref, alt), names_to="type", values_to="cDNA") %>%
        mutate(pep_id = ifelse(type == "alt", mut_id, sub("([0-9]+)[a-zA-Z]+$", "\\1", mut_id))) %>%
        pep_tile(...)

    if (nrow(fus) > 0) {
        fdf = as_tibble(as.data.frame(fus))
        ref1 = fdf %>% mutate(gene_name=sub("(.*)-.*", "\\1", fusion)) %>%
            select(fusion, gene_name, gene_id=gene_id_5p, tx_id=tx_id_5p, cDNA=ref_nuc_5p)
        ref2 = fdf %>% mutate(gene_name=sub(".*-(.*)", "\\1", fusion)) %>%
            select(fusion, gene_name, gene_id=gene_id_3p, tx_id=tx_id_3p, cDNA=ref_nuc_3p)
        refs = bind_rows(ref1, ref2) %>% mutate(type="ref", pep_id=sub("fs$", "", gene_name))
        alt = fdf %>%
            mutate(gene_id = paste(gene_id_5p, gene_id_3p, sep=";"),
                   tx_id = paste(tx_id_5p, tx_id_3p, sep=";"),
                   gene_name = fusion,
                   type = "alt",
                   pep_id = gene_name) %>% #TODO: should add position
            select(fusion, gene_id, tx_id, type, pep_id, cDNA=alt_nuc)
        ctx = bind_rows(refs, alt) %>%
            mutate(var_id=fusion, mut_id=fusion) %>%
            arrange(fusion, type, pep_id)
        pep = bind_rows(pep, pep_tile(ctx, ...))
    }

    list(`All Variants` = vars %>% select(-CDSLOC) %>% gr2df() %>%
            dplyr::select(-tx_name, -(ref_nuc:alt_prot)) %>% distinct(),
         `Unique Protein-Coding` = gr2df(subs) %>% select(var_id, mut_id, everything()),
         `RNA Fusions` = as.data.frame(fus),
         `93 nt Peptides` = pep %>% select(var_id, mut_id, pep_id,
            gene_id:cDNA, n_tiles, BbsI_replaced, tiled, nt, peptide))
}
