#' Aggregate fusion VCFs into a table
#'
#' @param vr    A VRanges object with RNA fusions from readVcfAsRanges
#' @param txdb  A transcription database, eg. AnnotationHub()[["AH100643"]]
#' @param asm   A Genome sequence package object, eg. ::BSgenome.Hsapiens.NCBI.GRCh38
#' @param min_reads  The minimum number of split read support for a fusion
#' @param min_pairs  The minimum number of linked read support for a fusion
#' @param min_tools  The minimum number of tools that identify a fusion
#' @return      A DataFrame objects with fusions
#'
#' @export
fusion = function(vr, txdb, asm, min_reads=NULL, min_pairs=NULL, min_tools=NULL) {
    # GT: ref allele/i alt allele (always './1'?)
    # DV: split reads
    # RV: discordant mates supporting translocation
    # FFPM: fusion fragments per million RNA fragments
    if (!is.null(min_reads))
        vr = vr[vr$DV >= min_reads]
    if (!is.null(min_pairs))
        vr = vr[vr$DV + vr$RV >= min_reads + min_pairs]
    if (!is.null(min_tools))
        vr = vr[unlist(vr$TOOL_HITS) >= min_tools]
    if (length(vr) == 0)
        return(DataFrame())

    # each possible combination of left and right transcripts from break
    res = tx_combine_break(txdb, vr)
    if (is.null(res))
        return(DataFrame())

    res$alt_nuc = get_coding_seq(asm, txdb,
        subseq(res$ref_nuc_5p, 1, res$break_cdsloc_5p),
        subseq(res$ref_nuc_3p, res$break_cdsloc_3p),
        include_stop = FALSE
    )

    is_fs = (res$break_cdsloc_5p %% 3 - (res$break_cdsloc_3p-1) %% 3) != 0
    res$fusion[is_fs] = paste0(res$fusion[is_fs], "fs")

    res[!duplicated(res$ref_nuc_5p) | !duplicated(res$ref_nuc_3p) |
        !duplicated(res$alt_nuc),]
}

#' Combine break info from each possible left and right side transcript
#'
#' @param txdb  A transcription database, eg. `AnnotationHub()[["AH100643"]]`
#' @param vr    A VRanges object with RNA fusions from `readVcfAsVRanges`
#' @return      A DataFrame with fusion coordinates and sequence information
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames seqlevelsStyle
#' @importFrom GenomicFeatures transcripts cdsBy
#' @keywords internal
tx_combine_break = function(txdb, vr) {
    flabs = paste(unlist(vr$GENEA), unlist(vr$GENEB), sep="-")
    g1 = GRanges(seqnames(vr), ranges(vr), sapply(vr$ORIENTATION, `[`, i=1))
    alt_loc = strsplit(gsub("\\[|\\]|N", "", alt(vr)), ":")
    g2 = GRanges(sapply(alt_loc, `[`, i=1), sapply(alt_loc, `[`, i=2),
                 sapply(vr$ORIENTATION, `[`, i=2))
    seqlevelsStyle(g2) = seqlevelsStyle(g1)

    tx = transcripts(txdb)
    coding_ranges = cdsBy(txdb)
    res = rep(list(list()), length(g1))
    for (i in seq_along(g1)) {
        left = tx_by_break(asm, txdb, tx, coding_ranges, g1[i], type="left")
        right = tx_by_break(asm, txdb, tx, coding_ranges, g2[i], type="right")
        if (nrow(left) == 0 || nrow(right) == 0)
            next
        colnames(left) = paste0(colnames(left), "_5p")
        colnames(right) = paste0(colnames(right), "_3p")
        idx = expand.grid(l=seq_len(nrow(left)), r=seq_len(nrow(right)))
        res[[i]] = cbind(fusion=flabs[i], split_reads=vr$DV[i], split_pairs=vr$RV[i],
                         FFPM=vr$FFPM[i], left[idx$l,], right[idx$r,])
    }
    do.call(rbind, res[sapply(res, length) > 0])
}

#' Get transcript incl coords left or right of break
#'
#' @param txdb  A transcription database, eg. AnnotationHub()[["AH100643"]]
#' @param tx    A list of transcripts obtained from `transcripts(txdb)`
#' @param coding_ranges  A list of exon coordinates by gene from `cdsBy(txdb)`
#' @param gr    GenomicRanges object of break location
#' @param type  Whether we want info for the 'left' or 'right' side of the break
#' @return      A DataFrame with sequence information
#'
#' @keywords internal
tx_by_break = function(asm, txdb, tx, coding_ranges, gr, type="left") {
    exo = exonsByOverlaps(txdb, gr)
    txo = subsetByOverlaps(coding_ranges, gr) # cdsByOverlaps can not take txdb
    if (type == "left") {
        exon_bound_is_break = ifelse(strand(exo) == "+", end(exo), start(exo)) == start(gr)
    } else if (type == "right") {
        exon_bound_is_break = ifelse(strand(exo) == "+", start(exo), end(exo)) == start(gr)
    }
    if (any(exon_bound_is_break))
        exo = exo[exon_bound_is_break]
    cdss = coding_ranges[intersect(names(coding_ranges), names(txo))]
    cdss = cdss[sapply(cdss, function(x) any(exo$exon_id %in% x$exon_id))]

    seqs = extractTranscriptSeqs(asm, cdss)
    prot = suppressWarnings(translate(seqs))
    has_start = subseq(prot,1,1) == "M"
    has_stop = subseq(IRanges::reverse(prot),1,1) == "*"
    n_stop = vcountPattern("*", prot)
    seqs = seqs[has_start & has_stop & n_stop==1] #todo: alt init codons
    if (length(seqs) == 0)
        return(DataFrame())
    locs = unlist(genomeToTranscript(gr, txdb))[names(seqs)]
    locs2 = transcriptToCds(locs, txdb)

    re = DataFrame(mcols(locs2)[c("seq_name", "seq_strand", "seq_start", "tx_id")],
                   gene_id=tx$gene_id[match(names(locs2), names(tx))],
                   ref_nuc=seqs, break_cdsloc=start(locs2))
    if (type == "left") {
        re = re[!duplicated(subseq(re$ref_nuc, 1, re$break_cdsloc)),]
    } else if (type == "right") {
        re = re[!duplicated(subseq(re$ref_nuc, re$break_cdsloc)),]
    }
    re
}
