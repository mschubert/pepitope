#' Subset the peptide context for gene fusions
#'
#' @param res  A DataFrame object from `fusions`
#' @param ctx_codons  How many flanking codons each to include in the context
#' @return     A DataFrame object of gene fusions
#' @export
subset_context_fusion = function(res, ctx_codons) {
    break_codon_start_5p = floor((res$break_cdsloc_5p-1)/3) * 3 + 1
    ref_starts_5p = break_codon_start_5p - ctx_codons*3
    ref_ends_5p = break_codon_start_5p + (ctx_codons+1)*3 - 1
    shift_5p = pmax(0, 1-ref_starts_5p) - pmax(0, ref_ends_5p-nchar(res$ref_nuc_5p))
    ref_nuc_5p = subseq(res$ref_nuc_5p, pmax(1, ref_starts_5p+shift_5p),
                        pmin(nchar(res$ref_nuc_5p)-3, ref_ends_5p+shift_5p))
    stopifnot(shift_5p %% 3 == 0)
    #stopifnot all peptides in original translation

    break_codon_start_3p = floor((res$break_cdsloc_3p-1)/3) * 3 + 1
    ref_starts_3p = break_codon_start_3p - ctx_codons*3
    ref_ends_3p = break_codon_start_3p + (ctx_codons+1)*3 - 1
    shift_3p = pmax(0, 1-ref_starts_3p) - pmax(0, ref_ends_3p-nchar(res$ref_nuc_3p))
    ref_nuc_3p = subseq(res$ref_nuc_3p, pmax(1, ref_starts_3p+shift_3p),
                        pmin(nchar(res$ref_nuc_3p)-3, ref_ends_3p+shift_3p))
    stopifnot(shift_3p %% 3 == 0)

    # fill results, annotate frameshifts and remove context dups
    concat = res$alt_nuc
    ctx_len = (ctx_codons*2 + 1) * 3
    is_fs = (res$break_cdsloc_5p %% 3 - (res$break_cdsloc_3p-1) %% 3) != 0
    end_3p = ref_starts_5p + ctx_len - 1
    res$alt_shift = pmin(0, nchar(concat)-end_3p) + pmax(0, shift_5p)
    end_3p[is_fs] = nchar(concat[is_fs])
    res$ref_nuc_5p = ref_nuc_5p
    res$ref_nuc_3p = ref_nuc_3p
    res$alt_nuc = subseq(concat, pmax(1, ref_starts_5p + res$alt_shift), end_3p)
    stopifnot(res$alt_shift %% 3 == 0)
    stopifnot(nchar(res$alt_nuc) %% 3 == 0)

    res[!duplicated(res$ref_nuc_5p) | !duplicated(res$ref_nuc_3p) |
        !duplicated(res$alt_nuc),]
}
