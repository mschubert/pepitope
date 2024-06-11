context("fusions")

ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
tx = suppressWarnings(transcripts(ens106))
cds = suppressWarnings(cdsBy(ens106))

# fusions from Mitelman database: EIF2AK2-STRN, NAIP-OCLN, AGO2-PTK2, DNM2-ILF3
chr = c("2", "5", "8", "19")
from = c(37141553, 70983775, 140635485, 10718403)
to = c(36905618, 69534694, 140879637, 10670516)
dir1 = c("-", "-", "-", "+")
dir2 = c("-", "+", "-", "+")
geneA = c("EIF2AK2", "NAIP", "AGO2", "DNM2")
geneB = c("STRN", "OCLN", "PTK2", "ILF3")

fus = function() {
    gr = GRanges(paste(chr, from, sep=":"))
    gr$ref = "N"
    gr$alt = sprintf(c("[%s:%i[N", "N]%s:%i]", "[%s:%i[N", "N[%s:%i["), chr, to)
    gr$ORIENTATION = do.call(IRanges::CharacterList, mapply(c, dir1, dir2, SIMPLIFY=FALSE))
    as(gr, "VRanges")
}

vr = fus()
ex = extract_fusion_ranges(vr)
left = lapply(ex$left, cds_by_break, txdb=ens106, cds=cds, type="left")
right = lapply(ex$right, cds_by_break, txdb=ens106, cds=cds, type="right")

test_that("ranges are extracted correctly", {
    v_left = unlist(ex$left)
    v_right = unlist(ex$right)
    expect_equal(as.character(seqnames(ex$left)), chr)
    expect_equal(as.character(seqnames(ex$right)), chr)
    expect_equal(IRanges::start(v_left), from)
    expect_equal(IRanges::start(v_right), to)
    expect_equal(as.character(GenomicRanges::strand(v_left)), dir1)
    expect_equal(as.character(GenomicRanges::strand(v_right)), dir2)
})

test_that("extracted transcripts all overlap break site", {
    tx1 = lapply(ex$left, transcriptsByOverlaps, x=ens106)
    expect_true(all(names(left[[1]]) %in% names(tx1[[1]])))
    expect_true(all(names(left[[2]]) %in% names(tx1[[2]])))
    expect_true(all(names(left[[3]]) %in% names(tx1[[3]])))
    expect_true(all(names(left[[4]]) %in% names(tx1[[4]])))

    tx2 = lapply(ex$right, transcriptsByOverlaps, x=ens106)
    expect_true(all(names(right[[1]]) %in% names(tx2[[1]])))
    expect_true(all(names(right[[2]]) %in% names(tx2[[2]])))
    expect_true(all(names(right[[3]]) %in% names(tx2[[3]])))
    expect_true(all(names(right[[4]]) %in% names(tx2[[4]])))
})

test_that("extracted transcripts belong to known genes", {
    tx1 = lapply(ex$left, transcriptsByOverlaps, x=ens106)
    gname1 = sapply(tx1, function(x) unique(sub("-[0-9]+$", "", na.omit(x$tx_external_name))))
    expect_equal(gname1, geneA)

    tx2 = lapply(ex$right, transcriptsByOverlaps, x=ens106)
    gname2 = sapply(tx2, function(x) unique(sub("-[0-9]+$", "", na.omit(x$tx_external_name))))
    expect_equal(gname2, geneB)
})

test_that("NAIP-OCLN fusion breakpoint is correctly assembled", {
    seq1 = add_seq_info(ex$left[[2]], left[[2]]['ENST00000194097'], asm, ens106, tx)
    seq2 = add_seq_info(ex$right[[2]], right[[2]]['ENST00000355237'], asm, ens106, tx)
    res = tx_combine_breaks(vr[2], list(seq1), list(seq2))
    expect_equal(nrow(res), 1)

    nuc_5p = subseq(res$ref_nuc_5p, 1, res$break_cdsloc_5p)
    nuc_3p = subseq(res$ref_nuc_3p, res$break_cdsloc_3p)
    alt_nuc = get_coding_seq(asm, ens106, nuc_5p, nuc_3p, include_stop=FALSE)

    read_support = "CACCCTGCCTTCCCTGGAATCTCTTGAAGTCTCAGGGACAATCCAGTCACAAGGT"
    expect_false(grepl(read_support, res$ref_nuc_5p))
    expect_false(grepl(read_support, res$ref_nuc_3p))
    expect_true(grepl(read_support, alt_nuc))
})

test_that("integrated fusion run", {
    res = annotate_fusions(vr, ens106, asm)
    cds_match = cds[names(cds) %in% c(res$tx_id_5p, res$tx_id_3p)]
    expect_true(all(tx[names(cds_match)]$tx_biotype == "protein_coding"))

    tx_seq = extractTranscriptSeqs(asm, cds_match)
    expect_true(all(tx_seq[res$tx_id_5p] == res$ref_nuc_5p))
    expect_true(all(tx_seq[res$tx_id_3p] == res$ref_nuc_3p))
})

# test cases:
# - multiple 3'UTR exons: ENST00000569022
