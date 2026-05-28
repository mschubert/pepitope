#' Logical vector whether sequence has START and STOP codon
#'
#' @param x  A DNAStringSet or AAStringSet object
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
