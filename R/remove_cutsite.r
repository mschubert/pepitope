#' Remove a Restriction Enzyme cut site but keep AA
#'
#' @param nuc   cDNA nucleotide string
#' @param site  Recognition site to be replaced (fwd+rev comp)
#' @param seed  Set random seed to select same changes on multiple runs
#' @return      cDNA with minimal changes to no longer contain the cut site
#'
#' @importFrom dplyr rowwise
#' @importFrom Biostrings subseq getGeneticCode DNAString translate replaceAt
#'      DNAStringSet reverseComplement vmatchPattern vcountPattern
#' @export
remove_cutsite = function(nuc, site, seed=NULL) {
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
