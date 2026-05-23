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
                       read_structures="7B12M+T", verbose=FALSE)
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
