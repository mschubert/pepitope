#' Make results report to save as xlsx sheets (full, filtered, peptides)
#'
#' @param vr       A VRanges object from `readVcfAsVRanges`
#' @param min_cov  Minimum number of reads to span the ALT allele
#' @param min_af   Minimum allele frequency of the ALT allele
#' @param pass     Whether to only include softFilterMatrix PASS
#'
#' @export
filter_variants = function(vr, min_cov=2, min_af=0.05, pass=TRUE) {
    if (!is.null(pass)) {
        passes = apply(softFilterMatrix(vr), 1, all)
        NAs = is.na(passes)
        if (all(NAs))
            stop("Can not apply filter [pass=TRUE] because all values are NA")
        vr = vr[!NAs & passes == pass]
    }

    if (!is.null(min_cov)) {
        NAs = is.na(altDepth(vr))
        if (all(NAs))
            stop("Can not apply filter [altDepth>=min_cov] because all values are NA")
        vr = vr[!NAs & altDepth(vr) >= min_cov]
    }

    if (!is.null(min_af)) {
        NAs = is.na(vr$AF)
        if (all(NAs))
            stop("Can not apply filter [AF>=min_af] because all values are NA")
        vr = vr[!NAs & vr$AF >= min_af]
    }

    vr
}
