context("read_structure")

test_that("rs_annotate finds barcode positions and orientation", {
    samples = data.frame(sample_id="sample1", patient="pat1", rep="1", origin="lib", barcode="GGG")
    constructs = list(pat1 = data.frame(gene_name="gene", mut_id="gene", pep_id="pep",
                                        tiled=TRUE, barcode="TTAACG"))
    fq = tempfile(fileext=".fq")
    writeLines(c(
        "@read1", "GGGNNNNNCGTTAA", "+", "IIIIIIIIIIIIII",
        "@read2", "GGGNNNNNTTAACG", "+", "IIIIIIIIIIIIII"
    ), fq)

    res = .rs_annotate(fq, samples, constructs, nrec=1)

    expect_equal(length(res$reads), 1L)
    expect_equal(res$counts[["B"]][1], 1L)
    expect_equal(res$counts[["^M"]][9], 1L)
    expect_equal(res$structure, "3B5S6^M")
})

test_that("read structure inference includes each barcode type meeting the threshold", {
    counts = data.frame(B=integer(14), `^B`=integer(14), M=integer(14),
                        `^M`=integer(14), check.names=FALSE)
    for (sample_op in c("B", "^B")) {
        for (construct_op in c("M", "^M")) {
            scores = counts
            scores[[sample_op]][1] = 20L
            scores[[construct_op]][9] = 20L
            expect_equal(.rs_format(scores, 3, 6, n_reads=100),
                         paste0("3", sample_op, "5S6", construct_op))
        }
    }

    counts$B[1] = 20L
    counts[["^M"]][9] = 20L
    for (weak_type in c("B", "^M")) {
        scores = counts
        scores[[weak_type]][] = 19L * (scores[[weak_type]] > 0)
        expect_equal(.rs_format(scores, 3, 6, n_reads=100),
                     if (weak_type == "B") "8S6^M" else "3B")
        expect_equal(.rs_format(scores, 3, 6, n_reads=100, min_frac=0.19), "3B5S6^M")
    }

    counts$B[1] = 10L
    counts[["^B"]][1] = 10L
    expect_equal(.rs_format(counts, 3, 6, n_reads=100), "8S6^M")
    expect_equal(.rs_format(counts, 3, 6, n_reads=50), "3B5S6^M")

    counts[["^B"]][1] = 0L
    counts$B[2] = 10L
    expect_equal(.rs_format(counts, 3, 6, n_reads=100), "8S6^M")
    counts[,] = 0L
    expect_error(.rs_format(counts, 3, 6, n_reads=100), "below min_frac")
})

test_that("partial structures preserve positions and orientations when parsed", {
    for (op in c("B", "^B", "M", "^M")) {
        counts = data.frame(B=integer(14), `^B`=integer(14), M=integer(14),
                            `^M`=integer(14), check.names=FALSE)
        counts[[op]][5] = 20L
        type = sub("^\\^", "", op)
        width = if (type == "B") 3L else 6L
        structure = .rs_format(counts, 3, 6, n_reads=100)
        expect_equal(structure, paste0("4S", width, op))
        parsed = .rs_parse(structure)
        present = if (type == "B") parsed$sample else parsed$construct
        absent = if (type == "B") parsed$construct else parsed$sample
        expect_equal(present$start, 5L)
        expect_equal(present$width, width)
        expect_equal(present$revcomp, startsWith(op, "^"))
        expect_equal(nrow(absent), 0L)
    }
    expect_error(.rs_parse("3B3B"), "at most one sample")
    expect_error(.rs_parse("6M6M"), "at most one construct")
    expect_error(.rs_parse("14S"), "must contain a sample")
})

test_that("overlapping barcode positions raise errors while adjacent positions are valid", {
    counts = data.frame(B=integer(14), `^B`=integer(14), M=integer(14),
                        `^M`=integer(14), check.names=FALSE)
    counts$B[1] = 20L
    counts$M[3] = 20L
    expect_error(.rs_format(counts, 3, 6, n_reads=100), "barcode positions overlap")
    counts$M[3] = 0L
    counts$M[4] = 20L
    expect_equal(.rs_format(counts, 3, 6, n_reads=100), "3B6M")

    counts[,] = 0L
    counts$M[1] = 20L
    counts[["^B"]][6] = 20L
    expect_error(.rs_format(counts, 3, 6, n_reads=100), "barcode positions overlap")
    counts[["^B"]][6] = 0L
    counts[["^B"]][7] = 20L
    expect_equal(.rs_format(counts, 3, 6, n_reads=100), "6M3^B")
})

test_that("read annotation passes min_frac and uses the actual number of inspected reads", {
    samples = data.frame(sample_id="sample1", patient="pat1", rep="1",
                         origin="lib", barcode="GGG")
    constructs = list(pat1=data.frame(gene_name="gene", mut_id="gene", pep_id="pep",
                                     tiled=TRUE, barcode="TTAACG"))
    fq = tempfile(fileext=".fq")
    on.exit(unlink(fq))
    reads = c("GGGNNNNNCGTTAA", rep(strrep("N", 14), 4))
    writeLines(as.vector(rbind(paste0("@r", seq_along(reads)), reads,
                               "+", strrep("I", 14))), fq)

    expect_equal(.rs_annotate(fq, samples, constructs, nrec=100)$structure, "3B5S6^M")
    expect_warning(ann <- .rs_annotate(fq, samples, constructs, nrec=100, min_frac=0.21),
                   "below min_frac")
    expect_equal(ann$structure, "unknown")
    expect_length(ann$reads, 5)
    expect_equal(ann$counts$B[1], 1L)
    expect_equal(ann$counts[["^M"]][9], 1L)
    for (structure in c("3B", "8S6^M"))
        expect_error(count_fastq(fq, samples, constructs, read_structure=structure, verbose=FALSE),
                     "Counting requires both")
})
