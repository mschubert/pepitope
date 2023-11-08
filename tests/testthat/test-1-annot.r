context("annotation")

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

ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
tx = transcripts(ens106, filter=SeqNameFilter("1"))
tx_coding = names(tx)[tx$tx_biotype=="protein_coding"]

vr = unlist(GRangesList(
    SYN = make_vr("1", c(5969221, 10647824, 18885617), "G", "A"),
    SNP = make_vr("1", c(1285573, 3763346, 43406323), c("C", "G", "T"), c("T", "A", "C")),
    STOP = make_vr("1", c(17119295, 39284153, 43312467), c("C", "A", "C"), c("T", "T", "A")),
    INS = make_vr("1", c(33010831, 21176231, 149509502), c("A", "T", "A"), c("AGGATGT", "TGCC", "AGAAGAC")),
    DEL = make_vr("1", c(1708855, 51840391, 201386873), c("TTCCTCCTCC", "ATCT", "CCCA"), c("T", "A", "C")),
    FS = make_vr("1", c(153810271, 33013400, 159062696), c("TG", "G", "GT"), c("T", "GATGTC", "G"))
))
ann = suppressWarnings(annotate_coding(vr, ens106, asm, tx_coding))

test_that("synonymous", {
    syn = ann[names(ann) == "SYN"]
    expect_equal(unique(syn$CONSEQUENCE), "synonymous")
    expect_true(all(syn$ref_prot == syn$alt_prot))
    expect_true(all(syn$ref_nuc != syn$alt_nuc))
    expect_equal(nchar(syn$ref_nuc), nchar(syn$alt_nuc))
    expect_equal(syn$alt_nuc, replaceDNA(syn$ref_nuc, syn$CDSLOC, syn$varAllele))
    expect_equal(syn$REFAA, syn$VARAA)
    expect_true(all(syn$silent_start == 1))
    expect_true(all(syn$silent_end == 0))

    ctx = subset_context(syn, ctx_codons=15)
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
    snp = ann[names(ann) == "SNP"]
    expect_true(all(snp$CONSEQUENCE == "nonsynonymous"))
    expect_true(all(snp$ref_prot != snp$alt_prot))
    expect_true(all(snp$ref_nuc != snp$alt_nuc))
    expect_equal(nchar(snp$ref_nuc), nchar(snp$alt_nuc))
    expect_equal(snp[strand(snp) == "+"]$alt, as.character(snp[strand(snp) == "+"]$varAllele))
    expect_equal(snp[strand(snp) == "-"]$alt, as.character(Biostrings::complement(snp[strand(snp) == "-"]$varAllele)))
    expect_true(all(snp$alt_nuc == replaceDNA(snp$ref_nuc, snp$CDSLOC, snp$varAllele)))
    expect_true(all(snp$REFAA != snp$VARAA))
    expect_true(all(snp$alt_prot == replaceProt(snp$ref_prot, snp$PROTEINLOC, snp$VARAA)))
    expect_true(all(snp$silent_start == 0))
    expect_true(all(snp$silent_end == 0))

    ctx = subset_context(snp, ctx_codons=15)
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
    stop = ann[names(ann) == "STOP"]
    expect_true(all(stop$CONSEQUENCE == "nonsense"))
    expect_true(all(as.character(stop$VARAA) == "*"))
    expect_true(all(nchar(stop$alt_prot) < nchar(stop$ref_prot)))

    ctx = subset_context(stop, ctx_codons=15)
    expect_true(all(nchar(ctx$ref_nuc) == 93))
    expect_true(all(nchar(ctx$alt_nuc) == 45))
    expect_true(all(nchar(ctx$ref_prot) == 31))
    expect_true(all(nchar(ctx$alt_prot) == 15))
    expect_equal(as.character(subseq(ctx$ref_prot, 1, 15)), as.character(ctx$alt_prot))
})

test_that("in-frame insertion", {
    ins = ann[names(ann) == "INS"]
    expect_true(all(ins$CONSEQUENCE == "insertion"))
    expect_true(all(nchar(ins$alt_prot) > nchar(ins$ref_prot)))

    ctx = subset_context(ins, ctx_codons=15)
    expect_true(all(nchar(ctx$alt_prot) > nchar(ctx$ref_prot)))
    expect_true(all(ctx$REFAA == "*" | ctx$silent_start > 0))
})

test_that("in-frame deletion", {
    del = ann[names(ann) == "DEL"]
    expect_true(all(del$CONSEQUENCE == "deletion"))
    expect_true(all(nchar(del$alt_prot) < nchar(del$ref_prot)))

    ctx = subset_context(del, ctx_codons=15)
    expect_equal(nchar(ctx$ref_nuc), 90+width(vr[ctx$QUERYID])-1)
    expect_true(all(nchar(ctx$alt_nuc) == 90))
})

test_that("frameshift", {
    fs = ann[names(ann) == "FS"]
    expect_true(all(fs$CONSEQUENCE == "frameshift"))

    ctx = subset_context(fs, ctx_codons=15)
    expect_true(all(nchar(ctx$ref_nuc) %% 3 == 0))
    expect_true(all(nchar(ctx$alt_nuc) %% 3 == 0))
    expect_true(all(nchar(ctx$alt_nuc) == 93))
    expect_true(all(nchar(ctx$alt_prot) == 31))
})
