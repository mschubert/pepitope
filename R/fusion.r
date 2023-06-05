#' Aggregate fusion VCFs into a table
#'
#' @param
#' @return
#'
#' @importFrom VariantAnnotation readVcfAsVRanges
#' @importFrom GenomicRanges GRanges
fusion = function(vr, txdb, asm, tx_coding) {
    f = "/DATA/m.schubert/nfcore-results/2023-05_CTh-mixed/rnafusion/megafusion/NSCLC_57.vcf"
    vr = readVcfAsVRanges(f, "GRCh38")

    # GT: ref allele/i alt allele (always './1'?)
    # DV: split reads
    # RV: discordant mates supporting translocation
    # FFPM: fusion fragments per million RNA fragments

    # get the pos+orientation for both ends
    # assemble possible variants same as in annotate_coding

    # filter by at least 2 read support + 1 pairs
    vr = vr[with(vr, DV >= 1 & RV >= 2 & unlist(TOOL_HITS) >= 2)]

    g1 = GRanges(seqnames(vr), ranges(vr), sapply(vr$ORIENTATION, `[`, i=1))
    alt_loc = strsplit(gsub("\\[|\\]|N", "", alt(vr)), ":")
    g2 = GRanges(sapply(alt_loc, `[`, i=1), sapply(alt_loc, `[`, i=2),
                 sapply(vr$ORIENTATION, `[`, i=2))

    ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
#    seqlevelsStyle(ens106) = "UCSC"
    is_prot = GeneBiotypeFilter("protein_coding")
    is_chr = SeqNameFilter(c(1:22,'X','Y'))
    tx = transcripts(ens106, filter=c(is_prot, is_chr))
    coding_ranges = cdsBy(ens106, filter=c(is_prot, is_chr))
    ex = exons(ens106, filter=c(is_prot, is_chr))

#    library(StructuralVariantAnnotation)
#    gr = breakpointRanges(readVcf(f), inferMissingBreakends=TRUE)
   asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
#    extractReferenceSequence(gr, asm, anchoredBases=45)

    get_pairs = function(left, right) {
        t1 = subsetByOverlaps(tx, left)
        e1 = subsetByOverlaps(ex, left)
        exon_end_is_break = ifelse(strand(e1) == "+", end(e1), start(e1)) == start(left)
        if (any(exon_end_is_break))
            e1 = e1[exon_end_is_break]
        tx_left = coding_ranges[intersect(names(coding_ranges), names(t1))]
        tx_left = tx_left[sapply(tx_left, function(x) any(e1$exon_id %in% x$exon_id))]

        t2 = subsetByOverlaps(tx, right)
        e2 = subsetByOverlaps(ex, right)
        exon_start_is_break = ifelse(strand(e2) == "+", start(e2), end(e2)) == start(right)
        if (any(exon_end_is_break))
            e2 = e2[exon_start_is_break]
        tx_right = coding_ranges[intersect(names(coding_ranges), names(t2))]
        tx_right = tx_right[sapply(tx_right, function(x) any(e2$exon_id %in% x$exon_id))]

        list(tx_left, tx_right)
    }
    tx_pairs = mapply(get_pairs, as(g1, "GRangesList"), as(g2, "GRangesList"), SIMPLIFY=FALSE)


    left1 = tx_pairs[[1]][[1]]
    left_coord = genomeToTranscript(g1[1], ens106)
    stopifnot(length(left_coord) == 1)
    left_break = as(left_coord[[1]][names(left1)], "DataFrame")
    left_break$break_cdsloc = start(left_break$X)
    left_break$ref_nuc = extractTranscriptSeqs(asm, left1)
    left_break$X = NULL

    right1 = tx_pairs[[1]][[2]]
    extractTranscriptSeqs(asm, right1)


}
