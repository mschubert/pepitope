# from: https://github.com/Bioconductor/VariantAnnotation/blob/devel/R/methods-predictCoding.R
# solves: https://github.com/Bioconductor/VariantAnnotation/pull/74
# solves: https://github.com/Bioconductor/VariantAnnotation/issues/86

#' @importFrom GenomicFeatures mapToTranscripts
#' @importFrom Biostrings AAStringSet DNAStringSetList GENETIC_CODE
#' @importFrom XVector subseq subseq<-
#' @importFrom S4Vectors %in% DataFrame elementNROWS
#' @importFrom MatrixGenerics rowRanges
#' @importFrom methods callGeneric as is
setGeneric("predictCoding",
    signature=c("query", "subject", "seqSource", "varAllele"),
    function(query, subject, seqSource, varAllele, ...)
        standardGeneric("predictCoding")
)

setMethod("predictCoding", c("IntegerRanges", "ANY", "ANY", "DNAStringSet"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    callGeneric(as(query, "GRanges"), subject, seqSource, varAllele, ...,
                ignore.strand=ignore.strand) 
})

setMethod("predictCoding", c("CollapsedVCF", "ANY", "ANY", "missing"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    rd <- rowRanges(query)
    alt <- alt(query) 
    if (is(alt, "CharacterList")) {
        alt <- .toDNAStringSetList(alt)
        if (sum(elementNROWS(alt)) == 0L) {
            stop("No nucleotide ALT values were detected.")
        }
    }
    rd <- rep(rowRanges(query), elementNROWS(alt))
    res <- callGeneric(rd, subject, seqSource, unlist(alt, use.names=FALSE), 
                ..., ignore.strand=ignore.strand)
    ## adjust QUERYID for expansion of rowRanges
    res$QUERYID <- rep(seq_len(length(alt)),
                       elementNROWS(alt))[res$QUERYID]
    res 
})

setMethod("predictCoding", c("ExpandedVCF", "ANY", "ANY", "missing"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    if (is(alt(query), "CharacterList")) {
      stop("alt(query) must be a DNAStringSet (not a CharacterList)")
    }
    callGeneric(rowRanges(query), subject, seqSource, alt(query), ..., 
                ignore.strand=ignore.strand) 
})

setMethod("predictCoding", c("GRanges", "ANY", "ANY", "DNAStringSet"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    .predictCoding(query, subject, seqSource, varAllele, ...,
                   ignore.strand=ignore.strand)
})

setMethod("predictCoding", c("VRanges", "ANY", "ANY", "missing"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    varAllele <- alt(query)
    query <- as(query, "GRanges")
    if (!is(varAllele, "DNAStringSet")) {
        tryCatch({
            varAllele <- DNAStringSet(varAllele)
        }, error=function(e) {
            stop(paste0("attempt to coerce 'alt' to DNAStringSet failed with ",
                 "error: ", conditionMessage(e)))
        })
    }
    .predictCoding(query, subject, seqSource, varAllele, ...,
                   ignore.strand=ignore.strand)
})

.predictCoding <-
    function(query, subject, seqSource, varAllele, ..., 
             cache=new.env(parent=emptyenv()), ignore.strand=FALSE)
{
    stopifnot(length(varAllele) == length(query))
    if (!any(seqlevels(query) %in% seqlevels(subject)))
        warning("none of seqlevels(query) match seqlevels(subject)")

    if (!exists(".__init__", cache, inherits=FALSE)) {
        cache[["cdsbytx"]] <- cdsBy(subject)
        cache[["txbygene"]] <- transcriptsBy(subject, "gene")
        cache[[".__init__"]] <- TRUE
    }

    map <- data.frame(geneid=rep(names(cache[["txbygene"]]), 
                          elementNROWS(cache[["txbygene"]])),
                      txid=mcols(unlist(cache[["txbygene"]], 
                          use.names=FALSE))[["tx_id"]],
                      stringsAsFactors=FALSE)

    txlocal <- .predictCodingGRangesList(query, cache[["cdsbytx"]], 
                   seqSource, varAllele, ..., ignore.strand=ignore.strand)
    txid <- mcols(txlocal)$TXID 
    mcols(txlocal)$GENEID <- map$geneid[match(txid, map$txid)]
    txlocal
}

.predictCodingGRangesList <- function(query, cdsbytx, seqSource, varAllele, 
                                      ..., genetic.code=GENETIC_CODE,
                                      if.fuzzy.codon="error", 
                                      ignore.strand=FALSE)
{
    if (ignore.strand)
        strand(query) <- "*"

    ## variant location in cds region
    mcols(query) <- append(mcols(query), DataFrame(varAllele=varAllele))
    txlocal <- .localCoordinates(query, cdsbytx, ignore.strand=FALSE, ...)

    ## reverse complement "-" strand
    valid <- rep(TRUE, length(txlocal))
    nstrand <- as.vector(strand(txlocal) == "-")
    if (any(nstrand)) {
        va <- mcols(txlocal)$varAllele
        va[nstrand] <- reverseComplement(va[valid & nstrand])
        mcols(txlocal)$varAllele <- va
    }

    ## frameshift
    refwidth <- width(txlocal)
    altallele <- mcols(txlocal)$varAllele
    fmshift <- abs(width(altallele) - refwidth) %% 3 != 0 
    if (any(fmshift))
        valid[fmshift] <- FALSE

    ## zero-width
    zwidth <- width(altallele) == 0
    if (any(zwidth)) {
        warning("records with missing 'varAllele' were ignored")
        valid[zwidth] <- FALSE 
        fmshift[zwidth] <- FALSE
    }

    ## reference codon sequences
    altpos <- (start(mcols(txlocal)$CDSLOC) - 1L) %% 3L + 1L
    refCodon <- varCodon <- .getRefCodons(txlocal, altpos, seqSource, cdsbytx)

    ## allowed characters that can't be translated
    ## "N", ".", "+" and "-"
    pattern <- "N|\\.|\\+|\\-"
    altCheck <- grepl(pattern, as.character(altallele, use.names=FALSE))
    refCheck <- grepl(pattern, as.character(refCodon, use.names=FALSE))
    noTrans <- rep(FALSE, length(txlocal)) 
    noTrans[altCheck | refCheck] <- TRUE
    valid[noTrans] <- FALSE
    if (any(altCheck))
        warning("'varAllele' values with 'N', '.', '+' or '-'",
                " were not translated")
    if (any(refCheck))
        warning("reference codons with 'N', '.', '+' or '-'",
                " were not translated")

    ## substitute and translate
    refAA <- varAA <- AAStringSet(rep("", length(txlocal))) 
    if (any(valid)) {
        ## 2 genetic.code versions
        alt.init.codons <- attr(genetic.code, "alt_init_codons")
        gc <- genetic.code
        gc.no.alt.init.codons <- genetic.code
        attr(gc.no.alt.init.codons, "alt_init_codons") <- character()
 
        ## ignore alt_init_codons 
        subseq(varCodon, altpos, width=refwidth) <- altallele
        refAA[valid] <- translate(refCodon[valid], 
                                  genetic.code=gc.no.alt.init.codons,
                                  if.fuzzy.codon=if.fuzzy.codon)
        varAA <- AAStringSet(rep("", length(txlocal))) 
        varAA[valid] <- translate(varCodon[valid], 
                                  genetic.code=gc.no.alt.init.codons, 
                                  if.fuzzy.codon=if.fuzzy.codon)

        ## respect alt_init_codons at the start of the CDS
        cds.start <- start(txlocal$CDSLOC) == 1L
        varAA.to.modify <- varCodon %in% alt.init.codons & cds.start
        refAA.to.modify <- refCodon %in% alt.init.codons & cds.start
        if (any(cds.start)) {
            varAA[varCodon %in% alt.init.codons & cds.start] <- "M"
            refAA[refCodon %in% alt.init.codons & cds.start] <- "M"
        }
    }

    ## results
    nonsynonymous <- as.character(refAA) != as.character(varAA) 
    consequence <- rep("synonymous", length(txlocal))
    consequence[nonsynonymous] <- "nonsynonymous" 
    consequence[fmshift] <- "frameshift"
    consequence[nonsynonymous & (as.character(varAA) %in% "*")] <- "nonsense" 
    consequence[zwidth | noTrans] <- "not translated" 
    consequence <- factor(consequence) 
 
    mcols(txlocal) <- append(mcols(txlocal), 
        DataFrame(GENEID=rep(NA_character_, nrow(mcols(txlocal))),
                  CONSEQUENCE=consequence, 
                  REFCODON=refCodon, 
                  VARCODON=varCodon, 
                  REFAA=refAA, VARAA=varAA))
    txlocal 
}

.getRefCodons <- function(txlocal, altpos, seqSource, cdsbytx)
{ 
    ## adjust codon end for 
    ## - width of the reference sequence
    ## - position of alt allele substitution in the codon
    cstart <- ((start(mcols(txlocal)$CDSLOC) - 1L) %/% 3L) * 3L + 1L
    cend <- cstart + (((altpos + width(txlocal) - 2L) %/% 3L) * 3L + 2L)
    txord <- match(mcols(txlocal)$TXID, names(cdsbytx))
    txseqs <- extractTranscriptSeqs(seqSource, cdsbytx[txord])
    DNAStringSet(substring(txseqs, cstart, cend))
}

.localCoordinates <- function(from, to, ignore.strand, ...)
{
    ## 'to' is a GRangesList of cds by transcript
    map <- mapToTranscripts(unname(from), to, ignore.strand=ignore.strand, ...)
    if (length(map) == 0) {
        res <- GRanges()
        mcols(res) <- DataFrame(REF=DNAStringSet(), ALT=DNAStringSetList(),
                                varAllele=DNAStringSet(), CDSLOC=IRanges(),
                                PROTEINLOC=IntegerList(), QUERYID=integer(),
                                TXID=character(), CDSID=IntegerList())
        return(res)
    }

    xHits <- map$xHits
    txHits <- map$transcriptsHits
    flat_to <- unlist(to) ## names needed for mapping

    ## FIXME: cdsid is expensive
    cdsid <- IntegerList(integer(0))
    map2 <- mapToTranscripts(unname(from)[xHits], flat_to,
                             ignore.strand=ignore.strand)
    cds <- mcols(flat_to)$cds_id[map2$transcriptsHits]
    ## CodingVariants() must fall within a coding region.
    ## mapToTranscripts() will map ranges that span intron 
    ## junctions and overlap multiple exons/cds regions. Because 
    ## of this, it's possible for 'map' to return a hit for a 
    ## query that 'map2' will not. (The subject in
    ## 'map' is a GRangesList and in 'map2' it's unlisted.)
    ## Only ranges identified by 'map' and 'map2' are kept.
    ## Ranges identified by 'map' only are discarded.
    if (length(cds)) {
        cdsid <- unique(splitAsList(cds, map2$xHits))
        keep <- logical(length(xHits))
        keep[unique(map2$xHits)] <- TRUE
        if (any(!keep)) {
            map <- map[keep]
            txHits <- map$transcriptsHits
            xHits <- map$xHits
        }
    }
    if (is.null(txid <- names(to)[txHits]))
        txid <- NA_integer_

    ## protein coordinates
    pends <- c(ceiling(start(map)/3), ceiling(end(map)/3))
    plocs <- unique(IntegerList(split(pends, rep(seq_len(length(pends)/2)), 2)))

    res <- from[xHits]
    strand(res) <- strand(map)
    mcols(res) <- append(values(res), 
        DataFrame(CDSLOC=ranges(map), 
                  PROTEINLOC=plocs, 
                  QUERYID=xHits, 
                  TXID=txid, CDSID=cdsid))
    res 
}
