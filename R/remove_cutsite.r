#' Remove a Restriction Enzyme cut site but keep AA
#'
#' @param nuc   cDNA nucleotide string
#' @param site  Recognition site to be replaced (fwd+rev comp)
#' @param seed  Set random seed to select same changes on multiple runs
#' @return      cDNA with minimal changes to no longer contain the cut site
remove_cutsite = function(nuc, site, seed=NULL) {
    set.seed(seed)
    revtrans = function(aa) {
        aa_split = strsplit(as.character(aa), "+")[[1]]
        tr = Biostrings::getGeneticCode()
        lapply(aa_split, function(x) names(tr)[tr == x]) %>%
            do.call(tidyr::crossing, .) %>%
            rowwise() %>%
            purrr::pmap_chr(paste0)
    }
    alt_nuc = function(nuc, match) {
        subs = subseq(Biostrings::DNAString(nuc), IRanges::start(match), IRanges::end(match))
        possib = revtrans(Biostrings::translate(subs, no.init.codon=TRUE))
        lapply(possib, Biostrings::replaceAt, x=nuc, at=match)
    }

    rc = function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
    m = c(vmatchPattern(site, nuc), vmatchPattern(rc(site), nuc)) %>% unlist()
    if (length(m) == 0)
        return(nuc)

    IRanges::start(m) = floor((IRanges::start(m)-1)/3) * 3 + 1
    IRanges::end(m) = ceiling((IRanges::end(m))/3) * 3

    nucs = nuc = Biostrings::DNAStringSet(nuc)
    for (i in seq_along(m))
        nucs = lapply(nucs, alt_nuc, match=m[i]) %>% do.call(c, .)
    nucs = Biostrings::DNAStringSet(nucs)

    valid = vcountPattern(site, nucs) + vcountPattern(rc(site), nucs) == 0
    nchange = mapply(function(i) stringdist::stringdist(subseq(nuc, m[i]), subseq(nucs, m[i])),
                     i=seq_along(m)) %>% as.matrix() %>% rowSums()
    min_valid = which(valid & nchange == min(nchange[valid]))
    as.character(nucs[sample(min_valid, 1)])
}
