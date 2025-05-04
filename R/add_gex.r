#' Add gene expression values to variant result
#'
#' @param res  Annotated variants from `annotate_coding()`
#' @param rec  Sample record from config file
#' @return     Annotated variants including gene counts and TPM
#'
#' @export
add_gex = function(res, rec) {
    counts = readr::read_tsv(rec$til_rna$count, col_names=c("gene_id", "gene_name", "count"))
    gex = readr::read_tsv(rec$til_rna$tpm) |>
        inner_join(counts) |>
        select(gene_id, gene_name, locus, count, TPM)

    res$rna_count = gex$count[match(res$GENEID, gex$gene_id)]
    res$rna_tpm = gex$TPM[match(res$GENEID, gex$gene_id)]
    res
}
