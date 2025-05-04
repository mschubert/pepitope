#' Remove a Restriction Enzyme cut site but keep AA in a tiled peptide data.frame
#'
#' @param pep   A data.frame of tiled peptides
#' @param ...   Named argumennts of cut sites, e.g. `BbsI="GAAGAC"`
#' @return      A data.frame with replace nucleotides and number of replacements
#' @export
remove_cutsite = function(pep, ...) {
    req = c("pep_id", "tiled")
    if (!all(req %in% colnames(pep)))
        stop("Required column(s) not found: ", paste(setdiff(req, colnames(pep)), collapse=", "))

    args = list(...)
    if (! length(args) == 1 && names(args) == "BbsI" && args[[1]] == "GAAGAC")
        stop("for now only BbsI implemented")

    pep$BbsI_replaced = vcountPattern("GAAGAC", pep$tiled) + vcountPattern("GTCTTC", pep$tiled)
    pep$tiled = sapply(pep$tiled, remove_cutsite_nuc, site="GAAGAC", seed=178529, USE.NAMES=FALSE)
#    stopifnot(pep$peptide == as.character(translate(DNAStringSet(pep$tiled), no.init.codon=TRUE)))
    pep
}

#' Remove a Restriction Enzyme cut site but keep AA
#'
#' @param nuc   cDNA nucleotide string
#' @param site  Recognition site to be replaced (fwd+rev comp)
#' @param seed  Set random seed to select same changes on multiple runs
#' @return      cDNA with minimal changes to no longer contain the cut site
#'
#' @importFrom Biostrings subseq getGeneticCode DNAString translate replaceAt
#'      DNAStringSet reverseComplement vmatchPattern vcountPattern
#' @keywords internal
remove_cutsite_nuc = function(nuc, site, seed=NULL) {
    set.seed(seed)
    revtrans = function(aa) {
        aa_split = strsplit(as.character(aa), "+")[[1]]
        tr = getGeneticCode()
        lapply(aa_split, function(x) names(tr)[tr == x]) |>
            do.call(tidyr::crossing, args=_) |>
            rowwise() |>
            purrr::pmap_chr(paste0)
    }
    alt_nuc = function(nuc, match) {
        subs = subseq(DNAString(nuc), IRanges::start(match), IRanges::end(match))
        possib = revtrans(translate(subs, no.init.codon=TRUE))
        lapply(possib, replaceAt, x=nuc, at=match)
    }

    rc = function(x) as.character(reverseComplement(DNAString(x)))
    m = c(vmatchPattern(site, nuc), vmatchPattern(rc(site), nuc)) |> unlist()
    if (length(m) == 0)
        return(nuc)

    IRanges::start(m) = floor((IRanges::start(m)-1)/3) * 3 + 1
    IRanges::end(m) = ceiling((IRanges::end(m))/3) * 3

    nucs = nuc = DNAStringSet(nuc)
    for (i in seq_along(m))
        nucs = lapply(nucs, alt_nuc, match=m[i]) |> do.call(c, args=_)
    nucs = DNAStringSet(nucs)

    valid = vcountPattern(site, nucs) + vcountPattern(rc(site), nucs) == 0
    nchange = mapply(function(i) stringdist::stringdist(subseq(nuc, m[i]), subseq(nucs, m[i])),
                     i=seq_along(m)) |> as.matrix() |> rowSums()
    min_valid = which(valid & nchange == min(nchange[valid]))
    as.character(nucs[sample(min_valid, 1)])
}
