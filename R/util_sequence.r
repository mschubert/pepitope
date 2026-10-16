#' Get the alternative coding sequence of a variant or gene fusion
#'
#' This takes into account any (newly introduced) STOP codons and UTR
#' readthroughs
#'
#' @param asm   Genomic sequence BSGenome object
#' @param txdb  A transcription database, eg. AnnotationHub()[["AH100643"]]
#' @param ...   A named DNAStringSet objects where each row is translated consecutively
#' @param include_stop  Whether to include the STOP codon
#' @return      A merged DNAStringSet object with the translated nucleotides
#'
#' @importFrom Biostrings xscat vmatchPattern translate
#' @importFrom GenomicFeatures threeUTRsByTranscript
#' @importFrom BSgenome getSeq
#' @keywords internal
get_coding_seq = function(asm, txdb, ..., include_stop=TRUE) {
    concat = xscat(...)
    if (length(concat) == 0)
        return(DNAStringSet())
    stops = suppressWarnings(vmatchPattern("*", translate(concat)))
    stop_offset = as.integer(!include_stop)
    stops = sapply(stops, function(s) (IRanges::start(s)[1]-stop_offset)*3)
    nostop = which(is.na(stops))

    # extend readthroughs to 3' UTR
    if (length(nostop) > 0) {
        utr3 = suppressWarnings(.txdb_cache_get(txdb, threeUTRsByTranscript(txdb)))
        tx_id_3p = names(rev(list(...))[[1]])
        if (length(intersect(names(utr3), tx_id_3p)) == 0)
            stop("3' DNAStringSet needs transcript IDs as names")

        for (i in nostop) {
            if (tx_id_3p[i] %in% names(utr3)) {
                exons_utr3 = getSeq(asm, utr3[[tx_id_3p[i]]])
                nuc_utr3 = do.call(xscat, as.list(exons_utr3))
                concat[i] = xscat(concat[i], nuc_utr3)
                pstop = suppressWarnings(vmatchPattern("*", translate(concat[i])))[[1]]
                stops[i] = (IRanges::start(pstop)[1]-stop_offset) * 3
            }
            if (is.na(stops[i]))
                stops[i] = floor(nchar(concat[i])/3) * 3
        }
    }

    stopifnot(stops %% 3 == 0)
    subseq(concat, 1, stops)
}

#' Logical vector whether sequence has START and STOP codon
#'
#' @param x  A DNAStringSet or AAStringSet object
#' @return  A logical vector indicating whether each sequence has a proper ORF
#'
#' @keywords internal
is_proper_orf = function(x) {
    if (inherits(x, "DNAStringSet"))
        x = suppressWarnings(translate(x))
    if (!inherits(x, "AAStringSet"))
        stop("'x' must be a DNAStringSet or AAStringSet object")

    has_start = subseq(x, 1, 1) == "M" #TODO: replace alt init codons by M?
    has_stop = subseq(IRanges::reverse(x), 1, 1) == "*"
    n_stop = vcountPattern("*", x)
    has_start & has_stop & n_stop == 1
}

#' @inheritParams is_proper_orf
#' @keywords internal
filter_proper_orf = function(x) {
    x[is_proper_orf(x)]
}
