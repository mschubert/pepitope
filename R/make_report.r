#' Make a variants report as named list of tables
#'
#' @param vars   Variant results
#' @param subs  Filtered variant results
#' @param fus   Fusion results
#' @param ...   Parameters passed to `pep_tile` (eg. tile_size, tile_ov)
#'
#' @importFrom dplyr mutate bind_rows desc distinct arrange as_tibble
#' @importFrom plyranges select
#' @export
make_report = function(vars, subs, fus=DataFrame(), ...) {
    gr2df = function(gr) as_tibble(as.data.frame(gr)) %>%
        select(var_id, everything()) %>%
        arrange(order(gtools::mixedorder(var_id)))

    # changes peptide, is unique and is expressed
    subs = subs[! subs$CONSEQUENCE %in% c("synonymous", "nonsense", "nostart")]
    alt_in_ref = function(a,r) grepl(as.character(a), as.character(r), fixed=TRUE)
    subs = subs[!mapply(alt_in_ref, a=subs$alt_prot, r=subs$ref_prot)]
    if ("rna_count" %in% colnames(S4Vectors::mcols(subs)) && !all(is.na(subs$rna_count)))
        subs = subs[!is.na(subs$rna_count) & subs$rna_count > 0] # & subs$rna_tpm > 0]
    subs = subs[!duplicated(subs$alt_prot)]

    # peptide is not contained within another peptide
    contained_in = function(i) any(grepl(ac2[i], ac2[-i], fixed=TRUE))
    ac = as.character(subs$alt_prot)
    ac2 = substr(ac, pmax(1, 1-subs$alt_shift/3), nchar(ac)-pmax(0, subs$alt_shift/3))
    any_con = sapply(seq_along(ac2), contained_in)
    subs = subs[!any_con]

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
            arrange(fusion, desc(type))
        pep = bind_rows(pep, pep_tile(ctx, ...))
    }

    list(`All Variants` = vars %>% select(-CDSLOC) %>% gr2df() %>%
            select(-tx_name, -(ref_nuc:alt_prot)) %>% distinct(),
         `Unique Protein-Coding` = gr2df(subs) %>% select(var_id, mut_id, everything()),
         `RNA Fusions` = as.data.frame(fus),
         `93 nt Peptides` = pep %>% select(var_id, mut_id, pep_id,
            gene_id:cDNA, n_tiles, BbsI_replaced, tiled, nt, peptide))
}
