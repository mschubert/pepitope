#' Aggregate fusion VCFs into a table
#'
#' @param vr    A VRanges object with RNA fusions from readVcfAsRanges
#' @param txdb  A transcription database, eg. AnnotationHub()[["AH100643"]]
#' @param min_reads  The minimum number of split read support for a fusion
#' @param min_pairs  The minimum number of linked read support for a fusion
#' @param min_tools  The minimum number of tools that identify a fusion
#' @param ctx_codons  Number of codonds for sequence context
#' @return      A DataFrame objects with fusions
#'
#' @importFrom VariantAnnotation readVcfAsVRanges
#' @importFrom GenomicRanges GRanges
#' @export
fusion = function(vr, txdb, asm, min_reads=NULL, min_pairs=NULL, min_tools=NULL, ctx_codons=15) {
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
    flabs = paste(unlist(vr$GENEA), unlist(vr$GENEB), sep="-")

    if (length(vr) == 0)
        return(DataFrame())

    g1 = GRanges(seqnames(vr), ranges(vr), sapply(vr$ORIENTATION, `[`, i=1))
    alt_loc = strsplit(gsub("\\[|\\]|N", "", alt(vr)), ":")
    g2 = GRanges(sapply(alt_loc, `[`, i=1), sapply(alt_loc, `[`, i=2),
                 sapply(vr$ORIENTATION, `[`, i=2))
    seqlevelsStyle(g2) = seqlevelsStyle(g1)

    tx = transcripts(txdb)
    coding_ranges = cdsBy(txdb)

    # from genomic break, get transcripts left and right + their coords
    tx_by_break = function(gr, type="left") {
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

    # each possible combination of left and right transcripts from break
    res = rep(list(list()), length(g1))
    for (i in seq_along(g1)) {
        left = tx_by_break(g1[i], type="left")
        right = tx_by_break(g2[i], type="right")
        if (nrow(left) == 0 || nrow(right) == 0)
            next
        colnames(left) = paste0(colnames(left), "_5p")
        colnames(right) = paste0(colnames(right), "_3p")
        idx = expand.grid(l=seq_len(nrow(left)), r=seq_len(nrow(right)))
        res[[i]] = cbind(fusion=flabs[i], split_reads=vr$DV[i], split_pairs=vr$RV[i],
                         FFPM=vr$FFPM[i], left[idx$l,], right[idx$r,])
    }
    res = do.call(rbind, res[sapply(res, length) > 0])
    if (is.null(res))
        return(DataFrame())

    # subset context
    break_codon_start_5p = floor((res$break_cdsloc_5p-1)/3) * 3 + 1
    ref_starts_5p = break_codon_start_5p - ctx_codons*3
    ref_ends_5p = break_codon_start_5p + (ctx_codons+1)*3 - 1
    shift_5p = pmax(0, 1-ref_starts_5p) - pmax(0, ref_ends_5p-nchar(res$ref_nuc_5p)+3)
    ref_nuc_5p = subseq(res$ref_nuc_5p, pmax(1, ref_starts_5p+shift_5p),
                        pmin(nchar(res$ref_nuc_5p), ref_ends_5p+shift_5p))
    stopifnot(shift_5p %% 3 == 0)
    #stopifnot all peptides in original translation

    break_codon_start_3p = floor((res$break_cdsloc_3p-1)/3) * 3 + 1
    ref_starts_3p = break_codon_start_3p - ctx_codons*3
    ref_ends_3p = break_codon_start_3p + (ctx_codons+1)*3 - 1
    shift_3p = pmax(0, 1-ref_starts_3p) - pmax(0, ref_ends_3p-nchar(res$ref_nuc_3p)+3)
    ref_nuc_3p = subseq(res$ref_nuc_3p, pmax(1, ref_starts_3p+shift_3p),
                        pmin(nchar(res$ref_nuc_3p), ref_ends_3p+shift_3p))
    stopifnot(shift_3p %% 3 == 0)

    ctx_len = (ctx_codons*2 + 1) * 3
    is_fs = (res$break_cdsloc_5p %% 3 - (res$break_cdsloc_3p-1) %% 3) != 0
    concat = xscat(subseq(res$ref_nuc_5p, 1, res$break_cdsloc_5p),
                   subseq(res$ref_nuc_3p, res$break_cdsloc_3p))
    end_3p = ref_starts_5p + ctx_len - 1
    bounded_end_3p = floor(pmin(nchar(concat), end_3p)/3) * 3
    stops = suppressWarnings(vmatchPattern("*", translate(concat)))
    stops = sapply(stops, function(s) (IRanges::start(s)[1]-1)*3)

    # extend transcripts w/o stops to UTRs
    nostop = which(is.na(stops))
    utr3 = threeUTRsByTranscript(txdb)
    for (i in nostop) {
        if (!res$tx_id_3p[i] %in% names(utr3))
            next
        nuc_utr3 = getSeq(asm, utr3[[res$tx_id_3p[i]]])
        concat[i] = xscat(concat[i], nuc_utr3)
        pstop = suppressWarnings(vmatchPattern("*", translate(concat[i])))[[1]]
        stops[i] = (IRanges::start(pstop)[1]-1) * 3
    }
    bounded_end_3p[is_fs] = stops[is_fs]

    # fill results, annotate frameshifts and remove context dups
    res$fusion[is_fs] = paste0(res$fusion[is_fs], "fs")
    res$ref_nuc_5p = ref_nuc_5p
    res$ref_nuc_3p = ref_nuc_3p
    res$alt_shift = pmin(0, bounded_end_3p-end_3p) + pmax(0, shift_5p)
    res$alt_nuc = subseq(concat, pmax(1, ref_starts_5p+res$alt_shift), bounded_end_3p)
    stopifnot(res$alt_shift %% 3 == 0)
    stopifnot(nchar(res$alt_nuc) %% 3 == 0)

    res[!duplicated(res$ref_nuc_5p) | !duplicated(res$ref_nuc_3p) |
        !duplicated(res$alt_nuc),]
}
