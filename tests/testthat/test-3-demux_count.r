context("demux_count")

example_barcodes = function(n, width=12) {
    alphabet = c("A", "C", "G", "T")
    vapply(seq_len(n) - 1L, function(i) {
        digits = (i %/% (4L ^ seq.int(0L, width - 1L))) %% 4L
        paste(alphabet[digits + 1L], collapse="")
    }, character(1))
}

valid_barcodes = example_barcodes(1000)
all_constructs = example_peptides(valid_barcodes)
sample_sheet = system.file("my_samples.tsv", package="pepitope")
samples = readr::read_tsv(sample_sheet)
fastq_file = example_fastq(sample_sheet, all_constructs, target_reads=100)

test_that("count source FASTQ directly", {
    dset = count_fastq(fastq_file, sample_sheet, all_constructs, valid_barcodes,
                       read_structure="7B12M+T", verbose=FALSE)
    expect_true(inherits(dset, "SummarizedExperiment"))
    expect_equal(samples$sample_id, colnames(dset))
    expect_equal(valid_barcodes, rownames(dset))

    n_reads = length(readLines(fastq_file)) / 4
    expect_equal(sum(SummarizedExperiment::assay(dset)), n_reads)
    expect_equal(SummarizedExperiment::colData(dset)$total_reads,
                 unname(colSums(SummarizedExperiment::assay(dset))))
    expect_equal(SummarizedExperiment::colData(dset)$mapped_reads,
                 unname(colSums(SummarizedExperiment::assay(dset))))
})

test_that("annotate_read_structure finds barcode positions and orientation", {
    samples = data.frame(sample_id="sample1", patient="pat1", rep="1", origin="lib", barcode="GGG")
    constructs = list(pat1 = data.frame(gene_name="gene", mut_id="gene", pep_id="pep",
                                        tiled=TRUE, barcode="TTAACG"))
    fq = tempfile(fileext=".fq")
    writeLines(c(
        "@read1", "GGGNNNNNCGTTAA", "+", "IIIIIIIIIIIIII",
        "@read2", "GGGNNNNNTTAACG", "+", "IIIIIIIIIIIIII"
    ), fq)

    res = annotate_read_structure(fq, samples, constructs, nrec=1)

    expect_equal(length(res$reads), 1L)
    expect_equal(res$counts[["B"]][1], 1L)
    expect_equal(res$counts[["M<"]][9], 1L)
    expect_equal(res$structure, "3B5S6M<")
})

test_that("count_fastq uses read structure orientation", {
    samples = data.frame(sample_id="sample1", patient="pat1", rep="1", origin="lib", barcode="GGG")
    constructs = list(pat1 = data.frame(gene_name="gene", mut_id="gene", pep_id="pep",
                                        tiled=TRUE, barcode="TTAACG"))
    fq = tempfile(fileext=".fq")
    writeLines(c("@read1", "GGGNNNNNCGTTAA", "+", "IIIIIIIIIIIIII"), fq)

    parsed = .rs_parse("3B5S6M<")
    expect_false(parsed$sample$revcomp)
    expect_true(parsed$construct$revcomp)

    dset = count_fastq(fq, samples, constructs, "TTAACG", read_structure="3B5S6M<", verbose=FALSE)
    expect_equal(SummarizedExperiment::assay(dset)[1, 1], 1)
})
