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
    flabs = paste(unlist(vr$GENEA), unlist(vr$GENEB), sep="-")

    g1 = GRanges(seqnames(vr), ranges(vr), sapply(vr$ORIENTATION, `[`, i=1))
    alt_loc = strsplit(gsub("\\[|\\]|N", "", alt(vr)), ":")
    g2 = GRanges(sapply(alt_loc, `[`, i=1), sapply(alt_loc, `[`, i=2),
                 sapply(vr$ORIENTATION, `[`, i=2))

    flt = ~ tx_biotype == "protein_coding" & SeqNameFilter(c(1:22,'X','Y'))
    ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
#    seqlevelsStyle(ens106) = "UCSC"
    tx = transcripts(ens106, filter=flt)
    coding_ranges = cdsBy(ens106, filter=flt)
    asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38

    # from genomic break, get transcripts left and right + their coords
    tx_by_break = function(gr, type="left") {
        exo = exonsByOverlaps(ens106, gr)
        txo = subsetByOverlaps(coding_ranges, gr) # cdsByOverlaps can not take ens106
        if (type == "left") {
            exon_bound_is_break = ifelse(strand(exo) == "+", end(exo), start(exo)) == start(gr)
        } else if (type == "right") {
            exon_bound_is_break = ifelse(strand(exo) == "+", start(exo), end(exo)) == start(gr)
        }
        if (any(exon_bound_is_break))
            exo = exo[exon_bound_is_break]
        cdss = coding_ranges[intersect(names(coding_ranges), names(txo))]
        cdss = cdss[sapply(cdss, function(x) any(exo$exon_id %in% x$exon_id))]
        if (length(cdss) == 0)
            return(DataFrame())

        seqs = extractTranscriptSeqs(asm, cdss)
        locs = unlist(genomeToTranscript(gr, ens106))[names(cdss)]
        locs2 = transcriptToCds(locs, ens106)

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
            return()
        colnames(left) = paste0(colnames(left), "_5p")
        colnames(right) = paste0(colnames(right), "_3p")
        idx = expand.grid(l=seq_len(nrow(left)), r=seq_len(nrow(right)))
        res[[i]] = cbind(fusion=flabs[i], left[idx$l,], right[idx$r,])
    }
    res = do.call(rbind, res[sapply(res, length) > 0])

    # subset context, annotate frameshifts, and remove context dups
}
