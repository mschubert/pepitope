#' Annotate VCF variants with coding changes
#'
#' @param vr    A VRanges object with SNVs and small indels
#' @param txdb  TxDb or EnsDb object
#' @param asm   Genomic sequence BSGenome object
#' @return      A GRanges object with annotated variants
#'
#' @examples
#' if (interactive()) {
#'     vcf = system.file("my_variants.vcf", package="pepitope")
#'     vr = readVcfAsVRanges(vcf)
#'     txdb = AnnotationHub::AnnotationHub()[["AH100643"]]
#'     asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
#'     annotate_coding(vr, txdb, asm)
#' }
#'
#' @importFrom Biostrings AAStringSet DNAStringSet DNAStringSetList GENETIC_CODE width
#' @importFrom Biostrings nchar replaceAt reverse reverseComplement subseq translate vcountPattern
#' @importFrom IRanges IRanges IntegerList
#' @importFrom S4Vectors %in% DataFrame elementNROWS mcols mcols<-
#' @importFrom GenomicFeatures cdsBy extractTranscriptSeqs genes mapToTranscripts transcriptsBy
#' @importFrom GenomicRanges strand strand<-
#' @importFrom VariantAnnotation sampleNames ref alt refDepth altDepth
#' @importFrom GenomeInfoDb seqnames seqnames<- genome genome<- isCircular isCircular<- seqinfo seqinfo<-
#' @importFrom methods as
#' @export
annotate_coding = function(vr, txdb, asm) {
    # workaround for EnsDb/BSgenome mismatch: https://github.com/jorainer/ensembldb/issues/159
    if ("UCSC" %in% seqlevelsStyle(vr) &&
            all(genome(txdb) == "GRCh38") && "UCSC" %in% seqlevelsStyle(txdb) &&
            all(genome(asm) == "hg38") && "UCSC" %in% seqlevelsStyle(asm)) {
        genome(seqinfo(asm))[] = "GRCh38"
    }
    if (length(intersect(seqlevels(vr), seqlevels(txdb))) == 0) {
        stop("No common `seqlevels()` in variants and transcript annotation. ",
            "Are you using the same Ensembl/UCSC styles?")
    }
    if (length(setdiff(seqlevels(vr), seqlevels(txdb))) > 0 &&
            length(setdiff(seqlevels(txdb), seqlevels(vr))) > 0) {
        stop("Incompatible `seqlevels()` in variants and transcript annotation. ",
            "Did you forget `chrs=\"default\"` in `filter_variants`?")
    }

    vr$sampleNames = sampleNames(vr)
    vr$AF = unlist(vr$AF)
    vr$ref = ref(vr)
    vr$alt = alt(vr)
    vr$cov_ref = refDepth(vr)
    vr$cov_alt = altDepth(vr)
    cds = .txdb_cache_get(txdb, cdsBy(txdb))
    codv = .coding_coordinates(as(vr, "GRanges"), cds)
    codv$varAllele = DNAStringSet(codv$alt)
    negative = strand(codv) == "-"
    codv$varAllele[negative] = reverseComplement(codv$varAllele[negative])

    tx_by_gene = .txdb_cache_get(txdb, transcriptsBy(txdb, "gene"))
    tx = unlist(tx_by_gene, use.names=FALSE)
    tx_match = match(codv$TXID, tx$tx_id)
    codv$GENEID = rep(names(tx_by_gene), elementNROWS(tx_by_gene))[tx_match]
    codv$tx_name = tx$tx_name[tx_match]
    gnames = .txdb_cache_get(txdb, genes(txdb))
    codv$gene_name = gnames$gene_name[match(codv$GENEID, gnames$gene_id)]

    # Extract each transcript once, then reuse its sequence for all overlapping variants.
    txids = unique(codv$TXID)
    ref_nuc = extractTranscriptSeqs(asm, cds[txids])
    ref_prot = suppressWarnings(translate(ref_nuc))
    tx_match = match(codv$TXID, txids)
    codv$ref_nuc = ref_nuc[tx_match]
    codv$ref_prot = ref_prot[tx_match]
    codv = codv[is_proper_orf(codv$ref_prot)]
    codv = .coding_changes(codv)
    codv$alt_nuc = get_coding_seq(asm, txdb, codv$alt_nuc)
    codv$alt_prot = suppressWarnings(translate(codv$alt_nuc))

    # name the variants
    codv$var_id = sprintf("%s:%i_%s/%s", seqnames(codv), IRanges::start(codv), codv$ref, codv$alt)
    pstarts = ceiling(start(codv$CDSLOC)/3) + codv$silent_start
    codv$mut_id = sprintf("%s_%s%i", codv$gene_name, codv$REFAA, pstarts)

    # check if we didn't change the length of any nuc
    stopifnot(with(codv,
        nchar(codv$ref_nuc) - nchar(codv$REFCODON) == nchar(codv$alt_nuc) - nchar(codv$VARCODON) |
        codv$CONSEQUENCE %in% c("frameshift", "nonsense", "nostart") |
        vcountPattern("*", codv$REFAA) > 0 | # transcript extension
        vcountPattern("*", codv$VARAA) > 0 # "missense"+nonsense
    ))

    infomap = match(seqnames(seqinfo(codv)), seqnames(seqinfo(gnames)))
    genome(codv) = genome(gnames)[infomap]
    isCircular(codv) = isCircular(gnames)[infomap]
    codv
}

# Coordinate mapping and local consequences adapted from VariantAnnotation's predictCoding:
# https://github.com/Bioconductor/VariantAnnotation/blob/devel/R/methods-predictCoding.R
# Includes the intron-spanning reference-width check discussed in:
# https://github.com/Bioconductor/VariantAnnotation/issues/86
.coding_coordinates = function(query, cds) {
    map = mapToTranscripts(unname(query), cds, ignore.strand=FALSE)
    if (length(map) == 0) {
        res = query[integer()]
        mcols(res) = cbind(mcols(res), DataFrame(CDSLOC=IRanges(), PROTEINLOC=IntegerList(),
            QUERYID=integer(), TXID=character()))
        return(res)
    }

    # Require a variant to map within an individual CDS, not merely across introns.
    flat_cds = unlist(cds[unique(map$transcriptsHits)])
    map_cds = mapToTranscripts(unname(query)[map$xHits], flat_cds, ignore.strand=FALSE)
    map = map[seq_along(map) %in% map_cds$xHits]

    protein_ends = c(ceiling(start(map)/3), ceiling(end(map)/3))
    protein_loc = unique(IntegerList(split(protein_ends, rep(seq_along(map), 2))))
    res = query[map$xHits]
    strand(res) = strand(map)
    mcols(res) = cbind(mcols(res), DataFrame(CDSLOC=ranges(map), PROTEINLOC=protein_loc,
        QUERYID=map$xHits, TXID=names(cds)[map$transcriptsHits]))
    res
}

.coding_changes = function(codv) {
    # Expand the affected reference interval to whole codons.
    codon_start = ((start(codv$CDSLOC) - 1L) %/% 3L) * 3L + 1L
    codon_end = ceiling(end(codv$CDSLOC)/3) * 3L
    codv$REFCODON = DNAStringSet(substring(codv$ref_nuc, codon_start, codon_end))
    invalid = width(codv) != width(codv$CDSLOC) | width(codv$REFCODON) < width(codv)
    if (any(invalid)) {
        warning("records with invalid reference width ignored")
        codv = codv[!invalid]
        codon_start = codon_start[!invalid]
        codon_end = codon_end[!invalid]
    }

    # Apply each allele once; derive its local codons from the same edited CDS.
    codv$alt_nuc = replaceAt(codv$ref_nuc, as(codv$CDSLOC, "IRangesList"),
        as(codv$varAllele, "DNAStringSetList"))
    delta = width(codv$varAllele) - width(codv)
    codv$VARCODON = DNAStringSet(substring(codv$alt_nuc, codon_start, codon_end + delta))
    frameshift = abs(delta) %% 3L != 0L
    missing = width(codv$varAllele) == 0L
    if (any(missing))
        warning("records with missing 'varAllele' were ignored")
    pattern = "N|\\.|\\+|\\-"
    alt_invalid = grepl(pattern, as.character(codv$varAllele))
    ref_invalid = grepl(pattern, as.character(codv$REFCODON))
    if (any(alt_invalid))
        warning("'varAllele' values with 'N', '.', '+' or '-' were not translated")
    if (any(ref_invalid))
        warning("reference codons with 'N', '.', '+' or '-' were not translated")
    valid = !frameshift & !missing & !alt_invalid & !ref_invalid

    codv$REFAA = codv$VARAA = AAStringSet(rep("", length(codv)))
    # A local codon window is internal to the protein unless it starts at the CDS start.
    codv$REFAA[valid] = translate(codv$REFCODON[valid], no.init.codon=TRUE)
    codv$VARAA[valid] = translate(codv$VARCODON[valid], no.init.codon=TRUE)
    at_start = valid & start(codv$CDSLOC) == 1L
    alt_init = attr(GENETIC_CODE, "alt_init_codons")
    codv$REFAA[at_start & codv$REFCODON %in% alt_init] = "M"
    codv$VARAA[at_start & codv$VARCODON %in% alt_init] = "M"

    nonsynonymous = codv$REFAA != codv$VARAA
    codv$CONSEQUENCE = rep("synonymous", length(codv))
    codv$CONSEQUENCE[nonsynonymous] = "nonsynonymous"
    codv$CONSEQUENCE[frameshift] = "frameshift"
    codv$CONSEQUENCE[nonsynonymous & codv$VARAA %in% "*"] = "nonsense"
    codv$CONSEQUENCE[missing | alt_invalid | ref_invalid] = "not translated"

    codv$silent_start = check_silent(codv$REFAA, codv$VARAA)
    codv$silent_end = check_silent(reverse(subseq(codv$REFAA, codv$silent_start+1)),
        reverse(subseq(codv$VARAA, codv$silent_start+1)))
    silent = codv$silent_start + codv$silent_end
    codv$CONSEQUENCE[start(codv$CDSLOC) == 1 & codv$VARCODON != "ATG"] = "nostart"
    codv$CONSEQUENCE[codv$REFAA == "*" & codv$VARAA != "*"] = "extension"
    codv$CONSEQUENCE[
        silent == nchar(codv$REFAA) & nchar(codv$VARAA) > nchar(codv$REFAA)] = "insertion"
    codv$CONSEQUENCE[
        silent == nchar(codv$VARAA) & nchar(codv$REFAA) > nchar(codv$VARAA)] = "deletion"
    codv
}

#' Adjust offsets for variants that start/end silently (no AA change)
#'
#' @param ref  Reference amino acid sequence
#' @param alt  Alternative amino acid sequence
#' @return     An integer vector of offsets
#'
#' @keywords internal
check_silent = function(ref, alt) {
    offsets = rep(0, length(ref))
    chk = seq_along(ref)
    for (i in seq_len(max(c(0, nchar(ref))))) {
        chk = chk[nchar(ref[chk]) >= i & nchar(alt[chk]) >= i]
        s_mtch = substr(ref[chk],i,i) == substr(alt[chk],i,i)
        if (! any(s_mtch))
            break
        chk = chk[s_mtch]
        offsets[chk] = offsets[chk] + 1
    }
    offsets
}
