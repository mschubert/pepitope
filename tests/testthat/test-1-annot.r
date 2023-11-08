context("annotation")

ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
tx = transcripts(ens106, filter=SeqNameFilter("1"))
tx_coding = names(tx)[tx$tx_biotype=="protein_coding"]

make_vr = function(chr, starts, ref, alt) {
    ends = starts + nchar(ref) - 1
    VariantAnnotation::VRanges(seqnames=chr, ranges=IRanges(starts, ends), ref=ref, alt=alt)
}

replaceDNA = function(what, where, to) {
    replaceAt(what, as(where, "IRangesList"), as(to, "DNAStringSetList"))
}

replaceProt = function(what, where, to) {
    replaceAt(what, as(unname(where), "IRangesList"), as(to, "AAStringSetList"))
}

test_that("synonymous", {
    vr = make_vr("1", c(5969221, 10647824, 18885617), "G", "A")

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    expect_equal(unique(ann$CONSEQUENCE), "synonymous")
    expect_true(all(ann$ref_prot == ann$alt_prot))
    expect_true(all(ann$ref_nuc != ann$alt_nuc))
    expect_equal(nchar(ann$ref_nuc), nchar(ann$alt_nuc))
    expect_equal(ann$alt_nuc, replaceDNA(ann$ref_nuc, ann$CDSLOC, ann$varAllele))
    expect_equal(ann$REFAA, ann$VARAA)
    expect_true(all(ann$silent_start == 1))
    expect_true(all(ann$silent_end == 0))

    ctx = subset_context(ann, ctx_codons=15)
    expect_true(all(ctx$ref_prot == ctx$alt_prot))
    expect_true(all(ctx$ref_nuc != ctx$alt_nuc))
    expect_true(all(nchar(ctx$ref_nuc) == 90))
    expect_true(all(nchar(ctx$alt_nuc) == 90))
    expect_equal(translate(ctx$ref_nuc), ctx$alt_prot)
    expect_equal(translate(ctx$alt_nuc), ctx$ref_prot)
    expect_true(all(nchar(ctx$ref_prot)), 30)
    expect_true(all(nchar(ctx$alt_prot)), 30)
})

test_that("SNPs", {
    vr = make_vr("1", c(1285573, 3763346, 43406323),
                 c("C", "G", "T"), c("T", "A", "C"))

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    expect_true(all(ann$CONSEQUENCE == "nonsynonymous"))
    expect_true(all(ann$ref_prot != ann$alt_prot))
    expect_true(all(ann$ref_nuc != ann$alt_nuc))
    expect_equal(nchar(ann$ref_nuc), nchar(ann$alt_nuc))
    expect_equal(ann[strand(ann) == "+"]$alt, as.character(ann[strand(ann) == "+"]$varAllele))
    expect_equal(ann[strand(ann) == "-"]$alt, as.character(Biostrings::complement(ann[strand(ann) == "-"]$varAllele)))
    expect_equal(ann$alt_nuc, replaceDNA(ann$ref_nuc, ann$CDSLOC, ann$varAllele))
    expect_true(all(ann$REFAA != ann$VARAA))
    expect_equal(ann$alt_prot, replaceProt(ann$ref_prot, ann$PROTEINLOC, ann$VARAA))
    expect_true(all(ann$silent_start == 0))
    expect_true(all(ann$silent_end == 0))

    ctx = subset_context(ann, ctx_codons=15)
    expect_true(all(ctx$ref_prot != ctx$alt_prot))
    expect_true(all(ctx$ref_nuc != ctx$alt_nuc))
    expect_true(all(nchar(ctx$ref_nuc) == 93))
    expect_true(all(nchar(ctx$alt_nuc) == 93))
    expect_equal(translate(ctx$ref_nuc), ctx$ref_prot)
    expect_equal(translate(ctx$alt_nuc), ctx$alt_prot)
    expect_true(all(nchar(ctx$ref_prot) == 31))
    expect_true(all(nchar(ctx$alt_prot) == 31))
})

test_that("nonsense", {
    vr = make_vr("1", c(17119295, 39284153, 43312467),
                 c("C", "A", "C"), c("T", "T", "A"))

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    expect_true(all(ann$CONSEQUENCE == "nonsense"))
    expect_true(all(as.character(ann$VARAA) == "*"))
    expect_true(all(nchar(ann$alt_prot) < nchar(ann$ref_prot)))

    ctx = subset_context(ann, ctx_codons=15)
    expect_true(all(nchar(ctx$ref_nuc) == 93))
    expect_true(all(nchar(ctx$alt_nuc) == 45))
    expect_true(all(nchar(ctx$ref_prot) == 31))
    expect_true(all(nchar(ctx$alt_prot) == 15))
    expect_equal(as.character(subseq(ctx$ref_prot, 1, 15)), as.character(ctx$alt_prot))
})

test_that("in-frame insertion", {
    vr = make_vr("1", c(33010831, 21176231, 149509502),
                 c("A", "T", "A"), c("AGGATGT", "TGCC", "AGAAGAC"))

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    expect_true(all(ann$CONSEQUENCE == "insertion"))
    expect_true(all(nchar(ann$alt_prot) > nchar(ann$ref_prot)))

    ctx = subset_context(ann, ctx_codons=15)
    expect_true(all(nchar(ctx$alt_prot) > nchar(ctx$ref_prot)))
    expect_true(all(ctx$REFAA == "*" | ctx$silent_start > 0))
})

test_that("in-frame deletion", {
    vr = make_vr("1", c(1708855, 51840391, 201386873),
                 c("TTCCTCCTCC", "ATCT", "CCCA"), c("T", "A", "C"))

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    expect_true(all(ann$CONSEQUENCE == "deletion"))
    expect_true(all(nchar(ann$alt_prot) < nchar(ann$ref_prot)))

    ctx = subset_context(ann, ctx_codons=15)
    expect_equal(nchar(ctx$ref_nuc), 90+width(vr[ctx$QUERYID])-1)
    expect_true(all(nchar(ctx$alt_nuc) == 90))
})

test_that("frameshift", {
    vr = make_vr("1", c(153810271, 33013400, 159062696),
                 c("TG", "G", "GT"), c("T", "GATGTC", "G"))

    ann = annotate_coding(vr, ens106, asm, tx_coding)
    expect_true(all(ann$CONSEQUENCE == "frameshift"))

    ctx = subset_context(ann, ctx_codons=15)
    expect_true(all(nchar(ctx$ref_nuc) %% 3 == 0))
    expect_true(all(nchar(ctx$alt_nuc) %% 3 == 0))
    expect_true(all(nchar(ctx$alt_nuc) == 93))
    expect_true(all(nchar(ctx$alt_prot) == 31))
})
