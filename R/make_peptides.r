#' Make a variants report as named list of tables
#'
#' @param subs  Variants within mutation context from `subset_context()`
#' @param fus   Variants within fusion context from `subset_context_fusions()`
#' @return      A data.frame with cDNA and peptide sequences
#'
#' @export
make_peptides = function(subs, fus=DataFrame()) {
    # changes peptide, is unique and is expressed
    subs = subs[! subs$CONSEQUENCE %in% c("synonymous", "nonsense", "nostart")]
    alt_in_ref = function(a,r) grepl(as.character(a), as.character(r), fixed=TRUE)
    subs = subs[!mapply(alt_in_ref, a=subs$alt_prot, r=subs$ref_prot)]
#    if ("rna_count" %in% colnames(S4Vectors::mcols(subs)) && !all(is.na(subs$rna_count)))
#        subs = subs[!is.na(subs$rna_count) & subs$rna_count > 0] # & subs$rna_tpm > 0]
    subs = subs[!duplicated(subs$alt_prot)]

    # peptide is not contained within another peptide
    contained_in = function(i) any(grepl(ac2[i], ac2[-i], fixed=TRUE))
    ac = as.character(subs$alt_prot)
    ac2 = substr(ac, pmax(1, 1-subs$alt_shift/3), nchar(ac)-pmax(0, subs$alt_shift/3))
    any_con = sapply(seq_along(ac2), contained_in)
    subs = subs[!any_con]

    pep = as.data.frame(subs) |>
        select(var_id, mut_id, gene_name, gene_id=GENEID, tx_id=tx_name,
               ref=ref_nuc, alt=alt_nuc) |>
        tidyr::pivot_longer(c(ref, alt), names_to="pep_type", values_to="cDNA") |>
        mutate(pep_id = ifelse(pep_type == "alt", mut_id, sub("([0-9]+)[a-zA-Z*]+$", "\\1", mut_id))) |>
        select(var_id, mut_id, pep_id, pep_type, everything())

    if (nrow(fus) == 0)
        return(pep)

    fdf = as_tibble(fus)
    ref1 = fdf |> mutate(gene_name=sub("(.*)--.*", "\\1", fusion)) |>
        select(fusion, gene_name, gene_id=gene_id_5p, tx_id=tx_id_5p, cDNA=ref_nuc_5p)
    ref2 = fdf |> mutate(gene_name=sub(".*--(.*)", "\\1", fusion)) |>
        select(fusion, gene_name, gene_id=gene_id_3p, tx_id=tx_id_3p, cDNA=ref_nuc_3p)
    refs = bind_rows(ref1, ref2) |> mutate(pep_type="ref", pep_id=sub("\\.fs$", "", gene_name))
    alt = fdf |>
        mutate(gene_id = paste(gene_id_5p, gene_id_3p, sep=";"),
               tx_id = paste(tx_id_5p, tx_id_3p, sep=";"),
               gene_name = fusion,
               pep_type = "alt",
               pep_id = gene_name) |> #TODO: should add position
        select(fusion, gene_id, tx_id, pep_type, pep_id, cDNA=alt_nuc)
    pep2 = bind_rows(refs, alt) |>
        mutate(var_id=fusion, mut_id=fusion) |>
        arrange(fusion, desc(pep_type)) |>
        select(-fusion)

    bind_rows(pep, pep2)
}
