#' Annotate VCF variants with coding changes
#'
#' @param vr    A VRanges object with SNVs and small indels
#' @param txdb  Txdb object
#' @param asm   Genomic sequence BSGenome object
#' @param tx_coding  Character vector of ENST0000 IDs that are protein coding
#' @param filter_variants  Apply soft filter matrix to the variants
#' @return      A GRanges object with annotated variants
#'
#' @importFrom GenomicFeatures transcripts threeUTRsByTranscript cdsBy
#'      extractTranscriptSeqs
#' @importFrom VariantAnnotation readVcfAsVRanges sampleNames ref alt refDepth
#'      altDepth predictCoding softFilterMatrix
#' @importFrom Biostrings subseq nchar reverse translate replaceAt DNAStringSet
#'      xscat vcountPattern
annotate_coding = function(vr, txdb, asm, tx_coding, filter_variants=FALSE) {
    if (filter_variants)
        vr = vr[apply(softFilterMatrix(vr), 1, all)]
    vr$sampleNames = sampleNames(vr)
    vr$AF = unlist(VR$AF)
    vr$ref = ref(vr)
    vr$alt = alt(vr)
    vr$cov_ref = refDepth(vr)
    vr$cov_alt = altDepth(vr)
    codv = predictCoding(vr, txdb, asm)

#    codv2 = predictCoding(vcf_tumor, txdb, asm) # 1637 var names @codv, 26k @codv2, 223 common?!
#    splice = locateVariants(vcf_diff, txdb, SpliceSiteVariants())

    tx = transcripts(txdb)
    utr3 = threeUTRsByTranscript(txdb)
    codv$tx_name = tx$tx_name[match(codv$TXID, tx$tx_id)]
    codv = codv[codv$tx_name %in% tx_coding]
    coding_ranges = cdsBy(txdb)[codv$TXID]
    codv$ref_nuc = unname(extractTranscriptSeqs(asm, coding_ranges))
    codv$ref_prot = translate(codv$ref_nuc)

    # filter for proper ORFs
    has_start = subseq(codv$ref_prot,1,1) == "M"
    has_stop = subseq(IRanges::reverse(codv$ref_prot),1,1) == "*"
    n_stop = vcountPattern("*", codv$ref_prot)
    codv = codv[has_start & has_stop & n_stop==1] #TODO: replace alt init codons by M

    # get coding sequences with updated variants
    codv$alt_nuc = replaceAt(codv$ref_nuc, as(codv$CDSLOC, "IRangesList"),
                             as(codv$varAllele, "DNAStringSetList"))

    # get protein sequences and adjust nuc for premature stop
    codv$alt_prot = translate(codv$alt_nuc)
    stops = vmatchPattern("*", codv$alt_prot)
    first = sapply(stops, function(s) IRanges::start(s)[1])
    changed = which(is.na(first) | first != nchar(codv$alt_prot))
    for (i in changed) {
        if (is.na(first[i])) {
            if (! codv$TXID[[i]] %in% names(utr3))
                next
            nuc_utr3 = getSeq(asm, utr3[[codv$TXID[i]]])
            codv$alt_nuc[i] = xscat(codv$alt_nuc[i], nuc_utr3)
            codv$alt_prot[i] = translate(codv$alt_nuc[i])
            first[i] = IRanges::start(vmatchPattern("*", codv$alt_prot[i])[[1]])[1]
        }
        codv$alt_prot[i] = subseq(codv$alt_prot[i], 1, first[i])
        codv$alt_nuc[i] = subseq(codv$alt_nuc[i], 1, first[i]*3)
    }

    check_silent = function(ref, alt) {
        offsets = rep(0, length(ref))
        chk = seq_along(ref)
        for (i in seq_len(max(nchar(ref)))) {
            chk = chk[nchar(ref[chk]) >= i & nchar(alt[chk]) >= i]
            s_mtch = substr(ref[chk],i,i) == substr(alt[chk],i,i)
            if (! any(s_mtch))
                break
            chk = chk[s_mtch]
            offsets[chk] = offsets[chk] + 1
        }
        offsets
    }
    codv$silent_start = check_silent(codv$REFAA, codv$VARAA)
    codv$silent_end = check_silent(reverse(subseq(codv$REFAA, codv$silent_start+1)),
                                   reverse(subseq(codv$VARAA, codv$silent_start+1)))
    silent = codv$silent_start + codv$silent_end

    codv$CONSEQUENCE = as.character(codv$CONSEQUENCE)
    codv$CONSEQUENCE[IRanges::start(codv$CDSLOC) == 1 & codv$VARCODON != "ATG"] = "nostart"
    codv$CONSEQUENCE[codv$REFAA == "*" & codv$VARAA != "*"] = "extension"
    codv$CONSEQUENCE[silent == nchar(codv$REFAA) & nchar(codv$VARAA) > nchar(codv$REFAA)] = "insertion"
    codv$CONSEQUENCE[silent == nchar(codv$VARAA) & nchar(codv$REFAA) > nchar(codv$VARAA)] = "deletion"

    codv$var_id = sprintf("%s:%i_%s/%s", GenomicRanges::seqnames(codv),
                          IRanges::start(codv), codv$ref, codv$alt)

    # check if we didn't change the length of any nuc
    stopifnot(with(codv,
        nchar(codv$ref_nuc) - nchar(codv$REFCODON) == nchar(codv$alt_nuc) - nchar(codv$VARCODON) |
        codv$CONSEQUENCE %in% c("frameshift", "nonsense", "nostart") |
        vcountPattern("*", codv$REFAA) > 0 | # transcript extension
        vcountPattern("*", codv$VARAA) > 0 # "missense"+nonsense
    ))
    codv
}
