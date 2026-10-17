test_that("sample sheets preserve identifiers, metadata names, and quoted text", {
    sheet = tempfile(fileext=".tsv")
    on.exit(unlink(sheet))
    writeLines(c(
        "sample_id\tminigenes\trep\torigin\tbarcode\textra field",
        "001\t007\t1\tLibrary\tAAC\t\"O'Brien #1\"",
        "\t\t\t\t\t",
        "002\t008\t2\tMock\tCCC\tNA"
    ), sheet)

    expect_warning(samples <- .read_samples(sheet), "Removing rows 2")
    expect_equal(samples$sample_id, c("001", "002"))
    expect_equal(samples$minigenes, c("007", "008"))
    expect_equal(samples$rep, c(1, 2))
    expect_equal(samples[["extra field"]], c("O'Brien #1", NA_character_))
})

test_that("barcode representation distinguishes construct and sample minigenes", {
    counts = matrix(seq_len(6), nrow=3,
                    dimnames=list(c("AAA", "CCC", "GGG"), c("s1", "s2")))
    bcs = data.frame(barcode=rownames(counts), minigenes=c("pat1", "common", "pat2"),
                      gene_name="GENE", mut_id="mut", pep_id="pep", pep_type="alt")
    meta = data.frame(sample_id=c("s1", "s2"), barcode=c("AAC", "TTG"),
                       minigenes=factor(c("pat1+common", "pat2")))
    reps = calc_representation(counts, bcs, meta)

    expect_equal(reps$minigenes, rep(bcs$minigenes, 2))
    expect_equal(as.character(reps$minigenes_sample), rep(c("pat1+common", "pat2"), each=3))
    expect_equal(reps$is_matched, c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))
    expect_equal(reps$value, as.vector(counts))
})

test_that("gene expression TSVs retain gene identifiers and numeric measurements", {
    counts = tempfile(fileext=".tsv")
    tpm = tempfile(fileext=".tsv")
    on.exit(unlink(c(counts, tpm)))
    writeLines(c("001\tGENE1\t1.5", "002\tGENE2\t2"), counts)
    writeLines(c("gene_id\tlocus\tTPM", "001\tchr1:1-2\t10", "002\tchr1:3-4\t20"), tpm)

    res = suppressMessages(add_gex(list(GENEID=c("002", "001", "003")),
                                   list(til_rna=list(count=counts, tpm=tpm))))
    expect_equal(res$rna_count, c(2, 1.5, NA_real_))
    expect_equal(res$rna_tpm, c(20, 10, NA_real_))
})

test_that("sample label templates use metadata columns literally", {
    meta = data.frame(sample_id=c("s1", "s2"), minigenes=factor(c("p1", "p2")),
                      rep=c(1, 2), origin=c("{minigenes}", "A\\B $1"))

    expect_equal(.format_sample_labels(meta, "{sample_id} ({minigenes})"), c("s1 (p1)", "s2 (p2)"))
    expect_equal(.format_sample_labels(meta, "{rep}: {origin} / {minigenes} / {minigenes}"),
                 c("1: {minigenes} / p1 / p1", "2: A\\B $1 / p2 / p2"))
    expect_equal(.format_sample_labels(meta, "All samples"), rep("All samples", 2))
    expect_equal(.format_sample_labels(meta[FALSE, ], "{sample_id}"), character())
    expect_error(.format_sample_labels(meta, "{missing}"), "Unknown sample label columns: missing")
})

test_that("coding sequence extension handles repeated transcripts, strands and missing UTRs", {
    asm = DNAStringSet(c(chr1="AACCCCCTAA", chr2="TTACCCCCC"))
    utr3 = GenomicRanges::GRangesList(
        forward = GRanges("chr1", IRanges(c(1, 8), width=3), strand="+"),
        reverse = GRanges("chr2", IRanges(c(7, 1), width=3), strand="-"),
        no_stop = GRanges("chr1", IRanges(1, width=3), strand="+")
    )
    testthat::local_mocked_bindings(.txdb_cache_get=function(txdb, value) utr3,
        .package="pepitope")
    left = DNAStringSet(setNames(rep("ATG", 6), rep("left", 6)))
    right = DNAStringSet(c(forward="CCC", forward="GGG", reverse="AAA",
        no_stop="CC", missing="CC", forward="TAA"))

    expect_equal(unname(as.character(get_coding_seq(asm, NULL, left, right))),
        c("ATGCCCAACTAA", "ATGGGGAACTAA", "ATGAAAGGGTAA", "ATGCCA", "ATG", "ATGTAA"))
    expect_equal(unname(as.character(get_coding_seq(asm, NULL, left, right, include_stop=FALSE))),
        c("ATGCCCAAC", "ATGGGGAAC", "ATGAAAGGG", "ATGCCA", "ATG", "ATG"))
    expect_equal(get_coding_seq(asm, NULL, DNAStringSet()), DNAStringSet())
    expect_error(get_coding_seq(asm, NULL, right[5]),
        "3' DNAStringSet needs transcript IDs as names", fixed=TRUE)
})
