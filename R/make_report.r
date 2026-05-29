#' Make a variants report as named list of tables
#'
#' @param vars  Variant results from `annotate_coding()`
#' @param subs  Variants within mutation context from `subset_context()`
#' @param fus   Variants within fusion context from `subset_context_fusions()`
#' @param tiled  A data.frame of the tiled peptide sequences
#' @return  A named list of data.frames for report output
#'
#' @examples
#' if (interactive()) {
#'     make_report(vars, subs, fus, tiled)
#' }
#'
#' @export
make_report = function(vars, subs, fus=DataFrame(), tiled) {
    gr2df = function(gr) {
        as_tibble(as.data.frame(gr)) |>
            select(var_id, everything()) |>
            arrange(order(gtools::mixedorder(var_id)))
    }

    mcols(vars) = mcols(vars)[names(mcols(vars)) != "CDSLOC"]

    list(`All Variants` = vars |> gr2df() |> select(-tx_name, -(ref_nuc:alt_prot)) |> distinct(),
         `Unique Protein-Coding` = gr2df(subs) |> select(var_id, mut_id, everything()),
         `RNA Fusions` = as.data.frame(fus),
         `93 nt Peptides` = tiled |> select(var_id, gene_name, mut_id, pep_id, pep_type,
            gene_id, tx_id, n_tiles, BbsI_replaced, tiled, nt, peptide))
}
