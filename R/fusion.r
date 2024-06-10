#' Aggregate fusion VCFs into a table
#'
#' @param vr    A VRanges object with RNA fusions from readVcfAsRanges
#' @param txdb  A transcription database, eg. AnnotationHub()[["AH100643"]]
#' @param asm   A Genome sequence package object, eg. ::BSgenome.Hsapiens.NCBI.GRCh38
#' @param min_reads  The minimum number of linked read support for a fusion
#' @param min_split_reads  The minimum number of split read support for a fusion
#' @param min_tools  The minimum number of tools that identify a fusion
#' @return      A DataFrame objects with fusions
#'
#' @importFrom GenomicFeatures transcripts cdsBy
#' @export
fusion = function(vr, txdb, asm, min_reads=NULL, min_split_reads=NULL, min_tools=NULL) {
    # GT: ref allele/i alt allele (always './1'?)
    # DV: split reads
    # RV: discordant mates supporting translocation
    # FFPM: fusion fragments per million RNA fragments
    if (!is.null(min_split_reads))
        vr = vr[vr$DV >= min_split_reads]
    if (!is.null(min_reads))
        vr = vr[vr$DV + vr$RV >= min_reads]
    if (!is.null(min_tools))
        vr = vr[unlist(vr$TOOL_HITS) >= min_tools]
    if (length(vr) == 0)
        return(DataFrame())

    # each possible combination of left and right transcripts from break
    tx = suppressWarnings(transcripts(txdb))
    cds = suppressWarnings(cdsBy(txdb))
    fr = extract_fusion_ranges(vr)
    left = lapply(fr$left, tx_by_break, asm=asm, txdb=txdb, tx=tx, cds=cds, type="left")
    right = lapply(fr$right, tx_by_break, asm=asm, txdb=txdb, tx=tx, cds=cds, type="right")
    res = tx_combine_breaks(vr, left, right)
    if (is.null(res))
        return(DataFrame())

    # extract fused coding sequence incl. possible STOP and 3'UTR extension
    nuc_5p = subseq(res$ref_nuc_5p, 1, res$break_cdsloc_5p)
    nuc_3p = subseq(res$ref_nuc_3p, res$break_cdsloc_3p)
    res$alt_nuc = get_coding_seq(asm, txdb, nuc_5p, nuc_3p, include_stop=FALSE)
    is_fs = (res$break_cdsloc_5p %% 3 - (res$break_cdsloc_3p-1) %% 3) != 0
    res$fusion[is_fs] = paste0(res$fusion[is_fs], "fs")

    res[!duplicated(res$ref_nuc_5p) | !duplicated(res$ref_nuc_3p) |
        !duplicated(res$alt_nuc),]
}

#' Convert a fusion VRanges object to left (5') and right (3') GRanges objects
#'
#' @param vr  A VRanges object with RNA fusions from `readVcfAsVRanges`
#' @return    A list of the 5' and 3' GRanges objects
#'
#' @importFrom S4Vectors splitAsList
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames seqlevelsStyle
#' @keywords internal
extract_fusion_ranges = function(vr) {
    g1 = GRanges(seqnames(vr), ranges(vr), sapply(vr$ORIENTATION, `[`, i=1))
    g2 = GRanges(gsub("\\[|\\]|N", "", alt(vr)), strand=sapply(vr$ORIENTATION, `[`, i=2))
    seqlevelsStyle(g2) = seqlevelsStyle(g1)
    list(left=splitAsList(g1), right=splitAsList(g2))
}

#' Combine break info from each possible left and right side transcript
#'
#' @param vr    A VRanges object with RNA fusions from `readVcfAsVRanges`
#' @param left  List of DataFrame objects containing the 5' of the fusion
#' @param right List of DataFrame objects containing the 3' of the fusion
#' @return      A DataFrame with fusion coordinates and sequence information
#'
#' @importFrom S4Vectors splitAsList
#' @keywords internal
tx_combine_breaks = function(vr, left, right) {
    combine_one = function(vr, left, right) {
        if (nrow(left) == 0 || nrow(right) == 0)
            return(DataFrame())
        lab = paste(vr$GENEA, vr$GENEB, sep="-")
        colnames(left) = paste0(colnames(left), "_5p")
        colnames(right) = paste0(colnames(right), "_3p")
        idx = expand.grid(l=seq_len(nrow(left)), r=seq_len(nrow(right)))
        cbind(fusion=lab, split_reads=vr$DV, split_pairs=vr$RV,
              FFPM=vr$FFPM, left[idx$l,], right[idx$r,])
    }
    res = mapply(combine_one, splitAsList(vr), left, right, SIMPLIFY=FALSE)
    do.call(rbind, res[sapply(res, length) > 0])
}

#' Get transcript incl coords left or right of break
#'
#' @param gr    GenomicRanges object of break location
#' @param asm   A Genome sequence package object, eg. ::BSgenome.Hsapiens.NCBI.GRCh38
#' @param txdb  A transcription database, eg. AnnotationHub()[["AH100643"]]
#' @param tx    A list of transcripts obtained from `transcripts(txdb)`
#' @param cds   A list of exon coordinates by gene from `cdsBy(txdb)`
#' @param type  Whether we want info for the 'left' or 'right' side of the break
#' @return      A DataFrame with sequence information
#'
#' @keywords internal
tx_by_break = function(gr, asm, txdb, tx, cds, type="left") {
    exo = exonsByOverlaps(txdb, gr)
    if (type == "left") {
        exon_bound_is_break = ifelse(strand(exo) == "+", end(exo), start(exo)) == start(gr)
    } else if (type == "right") {
        exon_bound_is_break = ifelse(strand(exo) == "+", start(exo), end(exo)) == start(gr)
    }
    if (any(exon_bound_is_break))
        exo = exo[exon_bound_is_break]

    txo = subsetByOverlaps(cds, gr) # cdsByOverlaps can not take txdb
    cdss = cds[intersect(names(cds), names(txo))]
    cdss = cdss[sapply(cdss, function(x) any(exo$exon_id %in% x$exon_id))]
    seqs = filter_proper_orf(extractTranscriptSeqs(asm, cdss))
    if (length(seqs) == 0)
        return(DataFrame())
    locs = unlist(genomeToTranscript(gr, txdb))[names(seqs)]
    locs2 = transcriptToCds(locs, txdb)

    DataFrame(mcols(locs2)[c("seq_name", "seq_strand", "seq_start", "tx_id")],
              gene_id=tx$gene_id[match(names(locs2), names(tx))],
              ref_nuc=seqs, break_cdsloc=start(locs2))
}
