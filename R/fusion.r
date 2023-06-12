#' Aggregate fusion VCFs into a table
#'
#' @param vr    A VRanges object with RNA fusions from readVcfAsRanges
#' @param txdb  A transcription database, eg. AnnotationHub()[["AH100643"]]
#' @param filter_fusions  Whether to only consider fusions with 2 tools, 2 reads
#' @return      A DataFrame objects with fusions
#'
#' @importFrom VariantAnnotation readVcfAsVRanges
#' @importFrom GenomicRanges GRanges
#' @export
fusion = function(vr, txdb, asm, filter_fusions=FALSE) {
    ens106 = txdb

    # GT: ref allele/i alt allele (always './1'?)
    # DV: split reads
    # RV: discordant mates supporting translocation
    # FFPM: fusion fragments per million RNA fragments
    if (filter_fusions)
        vr = vr[vr$DV >= 1 & vr$RV >= 2 & unlist(vr$TOOL_HITS) >= 2]
    flabs = paste(unlist(vr$GENEA), unlist(vr$GENEB), sep="-")

    g1 = GRanges(seqnames(vr), ranges(vr), sapply(vr$ORIENTATION, `[`, i=1))
    alt_loc = strsplit(gsub("\\[|\\]|N", "", alt(vr)), ":")
    g2 = GRanges(sapply(alt_loc, `[`, i=1), sapply(alt_loc, `[`, i=2),
                 sapply(vr$ORIENTATION, `[`, i=2))

    flt = ~ tx_biotype == "protein_coding" & SeqNameFilter(c(1:22,'X','Y'))
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

        seqs = extractTranscriptSeqs(asm, cdss)
        prot = suppressWarnings(translate(seqs))
        has_start = subseq(prot,1,1) == "M"
        has_stop = subseq(IRanges::reverse(prot),1,1) == "*"
        n_stop = vcountPattern("*", prot)
        seqs = seqs[has_start & has_stop & n_stop==1] #todo: alt init codons
        if (length(seqs) == 0)
            return(DataFrame())
        locs = unlist(genomeToTranscript(gr, ens106))[names(seqs)]
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
            next
        colnames(left) = paste0(colnames(left), "_5p")
        colnames(right) = paste0(colnames(right), "_3p")
        idx = expand.grid(l=seq_len(nrow(left)), r=seq_len(nrow(right)))
        res[[i]] = cbind(fusion=flabs[i], split_reads=vr$DV[i], split_pairs=vr$RV[i],
                         FFPM=vr$FFPM[i], left[idx$l,], right[idx$r,])
    }
    res = do.call(rbind, res[sapply(res, length) > 0])

    # subset context
    ctx_codons = 15
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

    is_fs = ((res$break_cdsloc_5p-1) %% 3 - (res$break_cdsloc_3p-1) %% 3) != 0
    concat = xscat(subseq(res$ref_nuc_5p, 1, res$break_cdsloc_5p),
                   subseq(res$ref_nuc_3p, res$break_cdsloc_3p)) # add UTRs?
    end_3p = pmin(nchar(concat), ref_starts_5p + (ctx_codons*2+1)*3)
    stops = suppressWarnings(vmatchPattern("*", translate(concat)))
    stops = sapply(stops, function(s) (IRanges::start(s)[1]-1)*3 + 1)
    end_3p[is_fs] = stops[is_fs]

    # fill results, annotate frameshifts and remove context dups
    res$fusion[is_fs] = paste0(res$fusion[is_fs], "fs")
    res$ref_nuc_5p = ref_nuc_5p
    res$ref_nuc_3p = ref_nuc_3p
    res$alt_shift = pmin(0, end_3p-ref_starts_5p - (ctx_codons*2+1)*3) + pmax(0, shift_5p)
    res$alt_nuc = subseq(concat, ref_starts_5p+res$alt_shift, end_3p)
    stopifnot(nchar(res$alt_nuc) %% 3 == 0)

    res[!duplicated(res$ref_nuc_5p) | !duplicated(res$ref_nuc_3p) |
        !duplicated(res$alt_nuc),]
}
