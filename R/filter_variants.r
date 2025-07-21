#' Make results report to save as xlsx sheets (full, filtered, peptides)
#'
#' @param vr       A VRanges object from `readVcfAsVRanges`
#' @param ...      Force filters by name (ignored)
#' @param min_cov  Minimum number of reads to span the ALT allele
#' @param min_af   Minimum allele frequency of the ALT allele
#' @param pass     Whether to only include softFilterMatrix PASS
#' @param sample   Only include if in `sampleNames(vr)` (required if more than one present)
#' @param chrs     Either "default" or a character vector of chromosome names
#'
#' @importFrom VariantAnnotation softFilterMatrix altDepth totalDepth
#' @importFrom Biobase sampleNames
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @export
filter_variants = function(vr, ..., min_cov=2, min_af=0.05, pass=TRUE, sample=NULL, chrs=NULL) {
    args = list(...)
    if (length(args) > 0)
        stop("Unknown filter arguments: ", paste(sQuote(args), collapse=", "))

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
        if (is.null(vr$AF)) {
            vec_af = altDepth(vr) / totalDepth(vr)
        } else if (inherits(vr$AF, "SimpleNumericList")) {
            if (any(sapply(vr$AF, length)) != 1)
                stop("'AF' field has multiple variants, convert to long format first")
            vec_af = unlist(vr$AF)
        } else {
            vec_af = vr$AF
        }

        NAs = is.na(vec_af)
        if (all(NAs))
            stop("Can not apply filter [AF>=min_af] because all values are NA")
        vr = vr[!NAs & vec_af >= min_af]
    }

    sn = sampleNames(vr)
    if (!is.null(sample)) {
        NAs = is.na(sn)
        if (all(NAs))
            stop("Can not apply filter [sampleNames] because all values are NA")
        vr = vr[!NAs & as.character(sn) %in% sample]
    } else if (length(unique(sn)) > 1) {
        stop("'vr' contains multiple samples but 'sample' argument not provided")
    }

    if (!is.null(chrs)) {
        if (length(chrs) == 1 && chrs == "default") {
            vr = keepStandardChromosomes(vr, pruning.mode="coarse")
        } else {
            vr = vr[seqnames(vr) %in% chrs]
        }
    }

    vr
}
