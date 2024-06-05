#' Annotate VCF variants with coding changes
#'
#' @param vr    A VRanges object with SNVs and small indels
#' @param txdb  TxDb or EnsDb object
#' @param asm   Genomic sequence BSGenome object
#' @return      A GRanges object with annotated variants
#'
#' @importFrom GenomicFeatures transcripts threeUTRsByTranscript cdsBy
#'      extractTranscriptSeqs
#' @importFrom VariantAnnotation readVcfAsVRanges sampleNames ref alt refDepth
#'      altDepth predictCoding
#' @importFrom Biostrings subseq nchar reverse translate replaceAt DNAStringSet
#'      xscat vcountPattern
#' @importFrom BSgenome getSeq
#' @importFrom GenomeInfoDb seqnames seqnames<- genome genome<-
#'      isCircular isCircular<- seqinfo seqinfo<-
#' @importFrom stringi stri_locate_first
#' @export
annotate_coding = function(vr, txdb, asm) {
    vr$sampleNames = sampleNames(vr)
    vr$AF = unlist(vr$AF)
    vr$ref = ref(vr)
    vr$alt = alt(vr)
    vr$cov_ref = refDepth(vr)
    vr$cov_alt = altDepth(vr)
    codv = predictCoding(vr, txdb, asm)

    gnames = genes(txdb)
    codv$gene_name = gnames$gene_name[match(codv$GENEID, gnames$gene_id)]

#    codv2 = predictCoding(vcf_tumor, txdb, asm) # 1637 var names @codv, 26k @codv2, 223 common?!
#    splice = locateVariants(vcf_diff, txdb, SpliceSiteVariants())

    tx = transcripts(txdb)
    codv$tx_name = tx$tx_name[match(codv$TXID, tx$tx_id)]
    if ("tx_biotype" %in% names(tx)) {
        tx_coding = names(tx)[tx$tx_biotype=="protein_coding"]
        codv = codv[codv$tx_name %in% tx_coding]
    }
    coding_ranges = cdsBy(txdb)[codv$TXID]
    codv$ref_nuc = extractTranscriptSeqs(asm, coding_ranges)
    codv$ref_prot = translate(codv$ref_nuc)

    # filter for proper ORFs
    has_start = subseq(codv$ref_prot,1,1) == "M"
    has_stop = subseq(IRanges::reverse(codv$ref_prot),1,1) == "*"
    n_stop = vcountPattern("*", codv$ref_prot)
    codv = codv[has_start & has_stop & n_stop==1] #TODO: replace alt init codons by M

    # get coding sequences with updated variants
    codv$alt_nuc = get_coding_seq(txdb,
        replaceAt(codv$ref_nuc, as(codv$CDSLOC, "IRangesList"),
                  as(codv$varAllele, "DNAStringSetList"))
    )
    codv$alt_prot = translate(codv$alt_nuc)

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

    # name the variants
    codv$var_id = sprintf("%s:%i_%s/%s", seqnames(codv), IRanges::start(codv), codv$ref, codv$alt)
    var_stop = stri_locate_first(codv$VARAA, fixed="*")[,1]
    vlabs = ifelse(is.na(var_stop), nchar(codv$VARAA), var_stop)
    mut_lab = ifelse(codv$CONSEQUENCE == "frameshift", "fs", substr(codv$VARAA, 1, vlabs))
    pstarts = unname(sapply(codv$PROTEINLOC, function(p) p[[1]])) + codv$silent_start
    codv$mut_id = sprintf("%s_%s%i%s", codv$gene_name, codv$REFAA, pstarts, mut_lab)

    # check if we didn't change the length of any nuc
    stopifnot(with(codv,
        nchar(codv$ref_nuc) - nchar(codv$REFCODON) == nchar(codv$alt_nuc) - nchar(codv$VARCODON) + 3 |
        codv$CONSEQUENCE %in% c("frameshift", "nonsense", "nostart") |
        vcountPattern("*", codv$REFAA) > 0 | # transcript extension
        vcountPattern("*", codv$VARAA) > 0 # "missense"+nonsense
    ))

    infomap = match(seqnames(seqinfo(codv)), seqnames(seqinfo(gnames)))
    genome(codv) = genome(gnames)[infomap]
    isCircular(codv) = isCircular(gnames)[infomap]
    codv
}
