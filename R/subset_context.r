#' Subset nucleotide/protein sequences to codon +/- 45 bp context
#'
#' @param codv  Annotated variants from `annotate_coding()`
#' @param ctx_codons  How many flanking codons each to include in the context
#' @return  GRanges object with sequence information of only context
#'
#' @importFrom dplyr `%>%`
#' @importFrom Biostrings vmatchPattern nchar
subset_context = function(codv, ctx_codons=15) {
    ctx = ctx_codons * 3
    stopAA = sapply(vmatchPattern("*", codv$VARAA), function(x) IRanges::start(x)[1]-1) * 3
    len_delta = pmin(nchar(codv$VARCODON), stopAA, na.rm=TRUE) - nchar(codv$REFCODON)

    roi_codon_start = floor((IRanges::start(codv$CDSLOC)-1)/3 + codv$silent_start) * 3 + 1
    roi_codon_end_ref = ceiling(IRanges::end(codv$CDSLOC)/3 - codv$silent_end) * 3
    len_ref = nchar(codv$ref_nuc) - 3
    len_alt = nchar(codv$alt_nuc) - 3

    is_frameshift = abs(len_delta) %% 3 != 0
    roi_codon_end_alt = pmin(ceiling((IRanges::end(codv$CDSLOC)+len_delta)/3 - codv$silent_end) * 3, len_alt)
    roi_codon_end_alt[is_frameshift] = nchar(codv$alt_nuc)[is_frameshift] - 3

    ctx_start = pmax(1, roi_codon_start - ctx)
    ctx_end_ref = pmin(len_ref, roi_codon_end_ref + ctx)
    ctx_end_alt = pmin(floor(len_alt/3)*3, roi_codon_end_alt + ctx) #TODO: better to force %%3 for nuc length?
    ctx_len_alt = ctx_end_alt - ctx_start + 1

    ctx_start_over = pmax(0, ctx + 1 - roi_codon_start)
    ctx_end_over = pmax(0, roi_codon_end_alt + ctx - len_alt)

    ext_by = pmin(len_alt, 2*ctx+3) - ctx_len_alt
    add_to_end = pmin(ext_by, ctx_start_over) %>% pmin(len_alt - ctx_start) %>% pmax(0)
    add_to_start = pmin(ext_by, ctx_end_over) %>% pmin(ctx_start-1) %>% pmax(0)
    add_to_start[codv$VARAA == "*"] = 0
    stopifnot(add_to_start == 0 | add_to_end == 0,
              codv$CONSEQUENCE == "extension" | (add_to_start %% 3 == 0 & add_to_end %% 3 == 0),
              c(add_to_start, add_to_end) <= ctx+3)

    # get nuc and protein sequences incl context
    codv$ref_nuc = subseq(codv$ref_nuc, ctx_start-add_to_start, ctx_end_ref+add_to_end) #TODO: merge if context overlaps? [only correct if same allele/read support]
    codv$alt_nuc = subseq(codv$alt_nuc, ctx_start-add_to_start, ctx_end_alt+add_to_end)
    codv$ref_prot = Biostrings::translate(codv$ref_nuc) #TODO: no alt init codons (atn M replaced by ATG above)
    codv$alt_prot = Biostrings::translate(codv$alt_nuc)
    codv$alt_nnuc = nchar(codv$alt_nuc)
    codv$alt_shift = add_to_end - add_to_start

    codv
}
