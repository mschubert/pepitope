import_package("GenomicFeatures", attach=TRUE)
import_package("VariantAnnotation", attach=TRUE)
import_package("dplyr", attach=TRUE)
sys = modules::import('sys')
subseq = Biostrings::subseq
nchar = Biostrings::nchar
reverse = Biostrings::reverse
vcountPattern = Biostrings::vcountPattern
vmatchPattern = Biostrings::vmatchPattern

#' Add gene expression values to variant result
#'
#' @param res  Annotated variants from `annotate_coding()`
#' @param rec  Sample record from config file
#' @return     Annotated variants including gene counts and TPM
add_gex = function(res, rec) {
    counts = readr::read_tsv(rec$til_rna$count, col_names=c("gene_id", "gene_name", "count"))
    gex = readr::read_tsv(rec$til_rna$tpm) %>%
        inner_join(counts) %>%
        dplyr::select(gene_id, gene_name, locus, count, TPM)

    res$rna_count = gex$count[match(res$GENEID, gex$gene_id)]
    res$rna_tpm = gex$TPM[match(res$GENEID, gex$gene_id)]
    res
}
