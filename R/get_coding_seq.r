#' Get the alternative coding sequence of a variant or gene fusion
#'
#' This takes into account any (newly introduced) STOP codons and UTR
#' readthroughs
#'
#' @param txdb  A transcription database, eg. AnnotationHub()[["AH100643"]]
#' @param ...   A named DNAStringSet objects where each row is translated consecutively
#' @return      A merged DNAStringSet object with the translated nucleotides
#'
#' @importFrom Biostrings xscat vmatchPattern translate
#' @importFrom GenomicFeatures threeUTRsByTranscript
#' @importFrom BSgenome getSeq
#' @keywords internal
get_coding_seq = function(txdb, ...) {
    concat = xscat(...)
    stops = suppressWarnings(vmatchPattern("*", translate(concat)))
    stops = sapply(stops, function(s) (IRanges::start(s)[1]-1)*3)
    nostop = which(is.na(stops))

    # extend readthroughs to 3' UTR
    utr3 = threeUTRsByTranscript(txdb)
    tx_id_3p = names(rev(list(...))[[1]])
    if (length(intersect(names(utr3), tx_id_3p)) == 0)
        stop("3' DNAStringSet needs transcript IDs as names")

    for (i in nostop) {
        if (tx_id_3p[i] %in% names(utr3)) {
            nuc_utr3 = getSeq(asm, utr3[[tx_id_3p[i]]])
            concat[i] = xscat(concat[i], nuc_utr3)
            pstop = suppressWarnings(vmatchPattern("*", translate(concat[i])))[[1]]
            stops[i] = (IRanges::start(pstop)[1]-1) * 3
        }
        if (is.na(stops[i]))
            stops[i] = floor(nchar(concat[i])/3) * 3
    }

    stopifnot(stops %% 3 == 0)
    subseq(concat, 1, stops)
}
