context("annotation")

ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
tx = transcripts(ens106, filter=SeqNameFilter("1"))
tx_coding = names(tx)[tx$tx_biotype=="protein_coding"]

make_vr = function(chr, starts, ref, alt) {
    ends = starts + nchar(ref) - 1
    VRanges(seqnames=chr, ranges=IRanges(starts, ends), ref=ref, alt=alt)
}

test_that("synonymous", {
    vr = make_vr("1", c(5969221, 10647824, 18885617), "G", "A")

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    with(ann, expect_equal(unique(CONSEQUENCE), "synonymous"))
    with(ann, expect_true(all(ref_prot == alt_prot)))
    with(ann, expect_true(all(ref_nuc != alt_nuc)))
    with(ann, expect_equal(nchar(ref_nuc), nchar(alt_nuc)))
    with(ann, expect_equal(alt_nuc,
        replaceAt(ref_nuc, as(CDSLOC, "IRangesList"), as(varAllele, "DNAStringSetList"))))
    with(ann, expect_equal(REFAA, VARAA))
    with(ann, expect_equal(unique(silent_start), 1))
    with(ann, expect_equal(unique(silent_end), 0))

    ctx = subset_context(ann, ctx_codons=15)
    with(ctx, expect_true(all(ref_prot == alt_prot)))
    with(ctx, expect_true(all(ref_nuc != alt_nuc)))
    with(ctx, expect_equal(unique(nchar(ref_nuc)), 90))
    with(ctx, expect_equal(unique(nchar(alt_nuc)), 90))
    with(ctx, expect_equal(translate(ref_nuc), alt_prot))
    with(ctx, expect_equal(translate(alt_nuc), ref_prot))
    with(ctx, expect_equal(unique(nchar(ref_prot)), 30))
    with(ctx, expect_equal(unique(nchar(alt_prot)), 30))
})

test_that("SNPs", {
    vr = make_vr("1", c(1285573, 3763346, 43406323),
                 c("C", "G", "T"), c("T", "A", "C"))

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    with(ann, expect_equal(unique(CONSEQUENCE), "nonsynonymous"))
    with(ann, expect_true(all(ref_prot != alt_prot)))
    with(ann, expect_true(all(ref_nuc != alt_nuc)))
    with(ann, expect_equal(nchar(ref_nuc), nchar(alt_nuc)))
    with(ann[strand(ann) == "+"], expect_equal(alt, as.character(varAllele)))
    with(ann[strand(ann) == "-"], expect_equal(alt, as.character(complement(varAllele))))
    with(ann, expect_equal(alt_nuc,
        replaceAt(ref_nuc, as(CDSLOC, "IRangesList"), as(varAllele, "DNAStringSetList"))))
    with(ann, expect_true(all(REFAA != VARAA)))
    with(ann, expect_equal(alt_prot,
        replaceAt(ref_prot, as(PROTEINLOC, "IRangesList"), as(VARAA, "AAStringSetList"))))
    with(ann, expect_equal(unique(silent_start), 0))
    with(ann, expect_equal(unique(silent_end), 0))

    ctx = subset_context(ann, ctx_codons=15)
    with(ctx, expect_true(all(ref_prot != alt_prot)))
    with(ctx, expect_true(all(ref_nuc != alt_nuc)))
    with(ctx, expect_equal(unique(nchar(ref_nuc)), 93))
    with(ctx, expect_equal(unique(nchar(alt_nuc)), 93))
    with(ctx, expect_equal(translate(ref_nuc), ref_prot))
    with(ctx, expect_equal(translate(alt_nuc), alt_prot))
    with(ctx, expect_equal(unique(nchar(ref_prot)), 31))
    with(ctx, expect_equal(unique(nchar(alt_prot)), 31))
})

test_that("nonsense", {
    vr = make_vr("1", c(17119295, 39284153, 43312467),
                 c("C", "A", "C"), c("T", "T", "A"))

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    with(ann, expect_equal(unique(CONSEQUENCE), "nonsense"))
    with(ann, expect_equal(unique(as.character(VARAA)), "*"))
    with(ann, expect_true(all(nchar(alt_prot) < nchar(ref_prot))))

    ctx = subset_context(ann, ctx_codons=15)
    with(ctx, expect_equal(unique(nchar(ref_nuc)), 93))
    with(ctx, expect_equal(unique(nchar(alt_nuc)), 45))
    with(ctx, expect_equal(unique(nchar(ref_prot)), 31))
    with(ctx, expect_equal(unique(nchar(alt_prot)), 15))
    with(ctx, expect_equal(as.character(subseq(ref_prot, 1, 15)), as.character(alt_prot)))
})

test_that("in-frame insertion", {
    vr = make_vr("1", c(33010831, 21176231, 149509502),
                 c("A", "T", "A"), c("AGGATGT", "TGCC", "AGAAGAC"))

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    with(ann, expect_equal(unique(CONSEQUENCE), "insertion"))
    with(ann, expect_true(all(nchar(alt_prot) > nchar(ref_prot))))

    ctx = subset_context(ann, ctx_codons=15)
    with(ctx, expect_true(all(nchar(alt_prot) > nchar(ref_prot))))
    with(ctx, expect_true(all(REFAA == "*" | silent_start > 0)))
})

test_that("in-frame deletion", {
    vr = make_vr("1", c(1708855, 51840391, 201386873),
                 c("TTCCTCCTCC", "ATCT", "CCCA"), c("T", "A", "C"))

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    with(ann, expect_equal(unique(CONSEQUENCE), "deletion"))
    with(ann, expect_true(all(nchar(alt_prot) < nchar(ref_prot))))

    ctx = subset_context(ann, ctx_codons=15)
    with(ctx, expect_equal(nchar(ref_nuc), 90+width(vr[QUERYID])-1))
    with(ctx, expect_equal(unique(nchar(alt_nuc)), 90))
})

test_that("frameshift", {
    vr = make_vr("1", c(153810271, 33013400, 159062696),
                 c("TG", "G", "GT"), c("T", "GATGTC", "G"))

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    with(ann, expect_equal(unique(CONSEQUENCE), "frameshift"))

    ctx = subset_context(ann, ctx_codons=15)
    with(ctx, expect_equal(unique(nchar(ref_nuc) %% 3), 0))
    with(ctx, expect_equal(unique(nchar(alt_nuc) %% 3), 0))
    with(ctx, expect_equal(unique(nchar(alt_nuc)), 93))
    with(ctx, expect_equal(unique(nchar(alt_prot)), 31))
})
