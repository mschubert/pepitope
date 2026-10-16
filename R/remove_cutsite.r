#' Remove restriction enzyme sites while preserving amino acids
#'
#' Minimizes nucleotide substitutions while preserving translation and excluding
#' all supplied recognition sites on both strands, including sites created by other edits.
#'
#' @param pep   A data.frame of tiled peptides
#' @param ...   Named recognition sequences containing A, C, G, and T, e.g.
#'      `BbsI="GAAGAC", EcoRI="GAATTC"`. All sites are removed together.
#' @return      A data.frame with synonymously recoded sequences and the number of
#'      original occurrences on either strand in one `<name>_replaced` column per site.
#'      Palindromic recognition sites are counted once per occurrence.
#'
#' @examples
#' pep = data.frame(pep_id="p1", tiled="ATGGAAGACGAATTC")
#' remove_cutsite(pep, BbsI="GAAGAC", EcoRI="GAATTC")
#'
#' @export
remove_cutsite = function(pep, ...) {
    req = c("pep_id", "tiled")
    if (!all(req %in% colnames(pep)))
        stop("Required column(s) not found: ", paste(setdiff(req, colnames(pep)), collapse=", "))

    args = list(...)
    if (length(args) == 0 || is.null(names(args)) || any(names(args) == "") ||
            anyDuplicated(names(args)))
        stop("Cut sites must have unique, nonempty names")
    sites = toupper(vapply(args, identity, character(1)))
    optimized = remove_cutsite_nuc(pep$tiled, site=sites, seed=178529)
    reverse = as.character(reverseComplement(DNAStringSet(sites)))

    for (i in seq_along(sites)) {
        patterns = unique(c(sites[i], reverse[i]))
        pep[[paste0(names(sites)[i], "_replaced")]] =
            Reduce(`+`, lapply(patterns, function(s) vcountPattern(s, pep$tiled)))
    }
    pep$tiled = optimized
    pep
}

#' Remove restriction enzyme sites while preserving amino acids
#'
#' Uses dynamic programming to minimize nucleotide substitutions while excluding
#' all recognition sites and their reverse complements everywhere in the result.
#' Incomplete terminal codons are kept unchanged.
#' Translation uses the standard genetic code with `no.init.codon=TRUE`.
#'
#' @param nuc   A character vector of cDNA nucleotide strings using unambiguous DNA bases
#' @param site  A character vector of recognition sequences to remove on both strands
#' @param seed  Seed for deterministic ordering of equally optimal synonymous codons
#' @return      cDNA sequences with the fewest nucleotide substitutions needed to
#'      remove the cut sites while preserving translation
#'
#' @importFrom Biostrings getGeneticCode DNAStringSet reverseComplement
#' @importFrom Biostrings subseq translate vcountPattern width
#' @keywords internal
remove_cutsite_nuc = function(nuc, site, seed=NULL) {
    if (!is.character(site) || length(site) == 0 || anyNA(site) ||
            any(!grepl("^[ACGT]+$", toupper(site))))
        stop("'site' must contain nonempty DNA recognition sequences using A, C, G, and T")
    site = toupper(site)
    sites = unique(c(site, as.character(reverseComplement(DNAStringSet(site)))))
    original = DNAStringSet(nuc)
    result = as.character(original)
    has_site = Reduce(`|`, lapply(sites, function(s) vcountPattern(s, original) > 0))
    if (!any(has_site))
        return(result)

    code = getGeneticCode()
    code = code[.pseudorandom_order(length(code), seed)]
    result[has_site] = remove_cutsite_nuc_cpp(result[has_site], sites, code)

    optimized = DNAStringSet(result[has_site])
    complete = 3L * (width(optimized) %/% 3L)
    before = translate(subseq(original[has_site], end=complete), no.init.codon=TRUE)
    after = translate(subseq(optimized, end=complete), no.init.codon=TRUE)
    stopifnot(identical(as.character(before), as.character(after)))
    for (s in sites)
        stopifnot(all(vcountPattern(s, optimized) == 0))
    result
}

.pseudorandom_order = function(n, seed=0) {
    if (n < 1)
        return(integer())

    seed = if (is.null(seed)) 0 else seed
    key = ((seq_len(n) + seed) * 1103515245 + 12345) %% 2^31
    order(key, seq_len(n))
}
