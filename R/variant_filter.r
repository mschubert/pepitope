#' Make results report to save as xlsx sheets (full, filtered, peptides)
#'
#' @param res      A full results GRanges object from `annotate_coding()`
#' @param min_cov  Minimum number of reads to span the ALT allele
#' @param min_af   Minimum allele frequency of the ALT allele
#'
#' @importFrom dplyr `%>%` rowwise mutate select arrange group_by ungroup
#'     tibble as_tibble n distinct
#' @importFrom Biostrings nchar translate DNAStringSet vcountPattern
#' @importFrom plyranges select
#' @export
variant_filter = function(res, min_cov=2, min_af=0.1) {
    # changes peptide, is unique and is expressed
    subs = subset_context(res[! res$CONSEQUENCE %in% c("synonymous", "nonsense", "nostart")])
    alt_in_ref = function(a,r) grepl(as.character(a), as.character(r), fixed=TRUE)
    subs = subs[!mapply(alt_in_ref, a=subs$alt_prot, r=subs$ref_prot)]
    if ("rna_count" %in% colnames(S4Vectors::mcols(subs)) && !all(is.na(subs$rna_count)))
        subs = subs[!is.na(subs$rna_count) & subs$rna_count > 0] # & subs$rna_tpm > 0]
    subs = subs[subs$AF >= min_af & subs$cov_alt >= min_cov]
    subs = subs[!duplicated(subs$alt_prot)]

    # peptide is not contained within another peptide
    contained_in = function(i) any(grepl(ac2[i], ac2[-i], fixed=TRUE))
    ac = as.character(subs$alt_prot)
    ac2 = substr(ac, pmax(1, 1-subs$alt_shift/3), nchar(ac)-pmax(0, subs$alt_shift/3))
    any_con = sapply(seq_along(ac2), contained_in)
    subs = subs[!any_con]

    # name the variants
    mut_lab = ifelse(subs$CONSEQUENCE == "frameshift", "fs", subs$VARAA)
    pstarts = unname(sapply(subs$PROTEINLOC, function(p) p[[1]])) + subs$silent_start
    subs$mut_id = sprintf("%s_%s%i%s", subs$gene_name, subs$REFAA, pstarts, mut_lab)

    table(nchar(subs$alt_nuc))
    stopifnot(nchar(subs$alt_nuc) %% 3 == 0)
    subs
}
