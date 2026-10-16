#' Add gene expression values to variant result
#'
#' @param res  Annotated variants from `annotate_coding()`
#' @param rec  Sample record from config file
#' @return     Annotated variants including gene counts and TPM
#'
#' @keywords internal
add_gex = function(res, rec) {
    counts = .read_tsv(rec$til_rna$count, header=FALSE,
                       col.names=c("gene_id", "gene_name", "count"),
                       colClasses=c(gene_id="character", gene_name="character"))
    gex = .read_tsv(rec$til_rna$tpm, colClasses=c(gene_id="character")) |>
        inner_join(counts) |>
        select(gene_id, gene_name, locus, count, TPM)

    res$rna_count = gex$count[match(res$GENEID, gex$gene_id)]
    res$rna_tpm = gex$TPM[match(res$GENEID, gex$gene_id)]
    res
}
