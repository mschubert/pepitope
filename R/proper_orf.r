#' Logical vector whether sequence has START and STOP codon
#'
#' @param x  A DNAStringSet or AAStringSet object
#'
#' @keywords internal
is_proper_orf = function(x) {
    UseMethod("is_proper_orf")
}

#' @inheritParams is_proper_orf
#' @keywords internal
is_proper_orf.DNAStringSet = function(x) {
    prot = suppressWarnings(translate(x))
    is_proper_orf(prot)
}

#' @inheritParams is_proper_orf
#' @keywords internal
is_proper_orf.AAStringSet = function(x) {
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
