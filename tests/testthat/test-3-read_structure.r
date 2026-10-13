context("read_structure")

test_that("rs_annotate finds barcode positions and orientation", {
    fq = tempfile(fileext=".fq")
    on.exit(unlink(fq))
    writeLines(c(
        "@read1", "GGGNNNNNCGTTAA", "+", "IIIIIIIIIIIIII",
        "@read2", "GGGNNNNNTTAACG", "+", "IIIIIIIIIIIIII"
    ), fq)

    res = .rs_annotate(fq, "GGG", "TTAACG", nrec=1)

    expect_equal(length(res$reads), 1L)
    expect_equal(res$counts[["B"]][1], 1L)
    expect_equal(res$counts[["^M"]][9], 1L)
    expect_equal(res$structure, "3B5S6^M")
})

test_that("read annotation counts forward and reverse barcode matches separately", {
    fq = tempfile(fileext=".fq")
    on.exit(unlink(fq))

    for (reverse in c(FALSE, TRUE)) {
        sample_seq = c(rep("AAC", 12), "GAG", "AAC")
        construct_seq = c(rep("TTAACG", 13), "CTGGAA")
        revcomp = c(rep(reverse, 8), rep(!reverse, 4), rep(reverse, 2))
        sample_seq[revcomp] = as.character(reverseComplement(DNAStringSet(sample_seq[revcomp])))
        construct_seq[revcomp] =
            as.character(reverseComplement(DNAStringSet(construct_seq[revcomp])))
        reads = paste0(sample_seq, "NNNNN", construct_seq)
        writeLines(as.vector(rbind(paste0("@r", seq_along(reads)), reads,
                                   "+", strrep("I", 14))), fq)

        counts = .rs_annotate(fq, "aac", "ttaacg", nrec=14)$counts
        main_ops = if (reverse) c("^B", "^M") else c("B", "M")
        other_ops = if (reverse) c("B", "M") else c("^B", "^M")
        for (i in seq_along(main_ops)) {
            position = c(1, 9)[i]
            expect_equal(counts[[main_ops[i]]][position], 9L)
            expect_equal(counts[[other_ops[i]]][position], 4L)
        }
    }
})

test_that("read matching extracts each barcode type in its inferred orientation", {
    for (sample_op in c("B", "^B")) {
        for (construct_op in c("M", "^M")) {
            sample_seq = DNAStringSet(c("AAC", "GAG", "AAC", "GAG"))
            construct_seq = DNAStringSet(c("TTAACG", "TTAACG", "CTGGAA", "CTGGAA"))
            if (sample_op == "^B")
                sample_seq = reverseComplement(sample_seq)
            if (construct_op == "^M")
                construct_seq = reverseComplement(construct_seq)
            reads = paste0(as.character(sample_seq), "NNNNN", as.character(construct_seq))
            reads = if (sample_op == "B") tolower(reads) else DNAStringSet(reads)

            matches = .rs_match(reads, paste0("3", sample_op, "5S6", construct_op),
                                 "aac", "ttaacg")
            sample = matches[matches$type == "B", ]
            construct = matches[matches$type == "M", ]
            expect_equal(sample$read, seq_len(4))
            expect_equal(construct$read, seq_len(4))
            expect_equal(sample$present, c(TRUE, FALSE, TRUE, FALSE))
            expect_equal(construct$present, c(TRUE, TRUE, FALSE, FALSE))
        }
    }
})

test_that("read matching uses the corresponding barcode list for each type", {
    matches = .rs_match(c("AACTTA", "TTAAAC"), "3B3M", "AAC", "TTA")
    for (type in c("B", "M"))
        expect_equal(matches$present[matches$type == type], c(TRUE, FALSE))
})

test_that("read matching includes only barcode types in partial structures", {
    cases = data.frame(read=c("NNGGGNNNNNNNNN", "NNCCCNNNNNNNNN",
                              "NNNNTTAACGNNNN", "NNNNCGTTAANNNN"),
                        structure=c("2S3B", "2S3^B", "4S6M", "4S6^M"),
                        type=c("B", "B", "M", "M"), start=c(3L, 3L, 5L, 5L))
    for (i in seq_len(nrow(cases))) {
        matches = .rs_match(c(cases$read[i], strrep("N", 14)), cases$structure[i],
                             "GGG", "TTAACG")
        expect_equal(nrow(matches), 2L)
        expect_equal(unique(matches$type), cases$type[i])
        expect_true(all(matches$start == cases$start[i]))
        expect_equal(matches$present, c(TRUE, FALSE))
    }
})

test_that("unknown structures produce no per-read barcode matches", {
    matches = .rs_match("GGGNNNNNTTAACG", "unknown", "GGG", "TTAACG")
    expect_equal(nrow(matches), 0L)
    expect_named(matches, c("read", "start", "width", "revcomp", "type", "present"))
    expect_type(matches$present, "logical")
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
    fq = tempfile(fileext=".fq")
    on.exit(unlink(fq))
    reads = c("GGGNNNNNCGTTAA", rep(strrep("N", 14), 4))
    writeLines(as.vector(rbind(paste0("@r", seq_along(reads)), reads,
                               "+", strrep("I", 14))), fq)

    expect_equal(.rs_annotate(fq, "GGG", "TTAACG", nrec=100)$structure, "3B5S6^M")
    expect_warning(ann <- .rs_annotate(fq, "GGG", "TTAACG", nrec=100, min_frac=0.21),
                   "below min_frac")
    expect_equal(ann$structure, "unknown")
    expect_length(ann$reads, 5)
    expect_equal(ann$counts$B[1], 1L)
    expect_equal(ann$counts[["^M"]][9], 1L)
})

test_that("read annotation warns and returns unknown for overlapping barcode positions", {
    fq = tempfile(fileext=".fq")
    on.exit(unlink(fq))
    writeLines(as.vector(rbind(paste0("@r", seq_len(5)), rep("GGGAAANNNNNNNN", 5),
                               "+", strrep("I", 14))), fq)

    expect_warning(ann <- .rs_annotate(fq, "GGG", "GGGAAA", nrec=5),
                   "barcode positions overlap")
    expect_equal(ann$structure, "unknown")
})
