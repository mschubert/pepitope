context("demux_count")

lib = "https://raw.githubusercontent.com/hawkjo/freebarcodes/master/barcodes/barcodes12-1.txt"
valid_barcodes = readr::read_tsv(lib, col_names=FALSE)$X1
all_constructs = example_peptides(valid_barcodes)
sample_sheet = system.file("my_samples.tsv", package="pepitope")
samples = readr::read_tsv(sample_sheet)
fastq_file = example_fastq(sample_sheet, all_constructs, target_reads=100)

test_that("demultiplex and count", {
    temp_dir = demux_fq(fastq_file, sample_sheet, read_structures="7B+T")
    fnames = list.files(temp_dir, pattern="\\.fq\\.gz$")
    expect_true(setequal(c(samples$sample_id, "unmatched"),
                         sub(".R1.fq.gz", "", fnames, fixed=TRUE)))

    dset = count_bc(temp_dir, all_constructs, valid_barcodes)
    expect_true(inherits(dset, "SummarizedExperiment"))
    expect_equal(samples$sample_id, colnames(dset))
    expect_equal(valid_barcodes, rownames(dset))
})
