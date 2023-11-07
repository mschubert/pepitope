context("annotation")

ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
tx = transcripts(ens106, filter=SeqNameFilter("1"))
tx_coding = names(tx)[tx$tx_biotype=="protein_coding"]

make_vr = function(chr, starts, ref, alt) {
    ends = starts + nchar(alt) - nchar(ref)
    VRanges(seqnames=chr, ranges=IRanges(starts, ends), ref=ref, alt=alt)
}

test_that("synonymous", {
    vr = make_vr("1", c(5969221, 10647824, 18885617), "G", "A")

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    with(ann, expect_true(unique(ref_prot == alt_prot)))
    with(ann, expect_true(unique(ref_nuc != alt_nuc)))
    with(ann, expect_equal(nchar(ref_nuc), nchar(alt_nuc)))
    with(ann, expect_equal(REFAA, VARAA))
    with(ann, expect_equal(unique(silent_start), 1))
    with(ann, expect_equal(unique(silent_end), 0))

    ctx = subset_context(ann, ctx_codons=15)
    with(ctx, expect_true(unique(ref_prot == alt_prot)))
    with(ctx, expect_true(unique(ref_nuc != alt_nuc)))
    with(ctx, expect_equal(unique(nchar(ref_nuc)), 90))
    with(ctx, expect_equal(unique(nchar(alt_nuc)), 90))
})

test_that("SNPs", {
    vr = make_vr("1", c(1285573, 3763346, 43406323),
                 c("C", "G", "T"), c("T", "A", "C"))

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    with(ann, expect_true(unique(ref_prot != alt_prot)))
    with(ann, expect_true(unique(ref_nuc != alt_nuc)))
    with(ann, expect_equal(nchar(ref_nuc), nchar(alt_nuc)))
    with(ann[strand(ann) == "+"], expect_equal(alt, as.character(varAllele)))
    with(ann[strand(ann) == "-"], expect_equal(alt, as.character(complement(varAllele))))
    with(ann, expect_equal(alt_nuc,
        replaceAt(ref_nuc, as(CDSLOC, "IRangesList"), as(varAllele, "DNAStringSetList"))))
    with(ann, expect_true(unique(REFAA != VARAA)))
    with(ann, expect_equal(alt_prot,
        replaceAt(ref_prot, as(PROTEINLOC, "IRangesList"), as(VARAA, "AAStringSetList"))))
    with(ann, expect_equal(unique(silent_start), 0))
    with(ann, expect_equal(unique(silent_end), 0))

    ctx = subset_context(ann, ctx_codons=15)
    with(ctx, expect_true(unique(ref_prot != alt_prot)))
    with(ctx, expect_true(unique(ref_nuc != alt_nuc)))
    with(ctx, expect_equal(unique(nchar(ref_nuc)), 93))
    with(ctx, expect_equal(unique(nchar(alt_nuc)), 93))
    with(ctx, expect_equal(translate(ref_nuc), ref_prot))
    with(ctx, expect_equal(translate(alt_nuc), alt_prot))
    with(ctx, expect_equal(unique(nchar(ref_prot)), 31))
    with(ctx, expect_equal(unique(nchar(alt_prot)), 31))
})
