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
        vr = vr[!is.na(passes) & passes == pass]
    }

    if (!is.null(min_cov))
        vr = vr[!is.na(altDepth(vr)) & altDepth(vr) >= min_cov]

    if (!is.null(min_af))
        vr = vr[!is.na(vr$AF) & vr$AF >= min_af]

    vr
}
