context("count_fastq")

valid_barcodes = example_barcodes(1000)
all_constructs = example_peptides(valid_barcodes)
sample_sheet = system.file("my_samples.tsv", package="pepitope")
samples = .read_samples(sample_sheet)
fastq_file = example_fastq(sample_sheet, all_constructs, target_reads=100)

test_that("count source FASTQ directly", {
    dset = count_fastq(fastq_file, sample_sheet, all_constructs, valid_barcodes,
                       read_structure="7B12M+T", verbose=FALSE)
    expect_true(inherits(dset, "SummarizedExperiment"))
    expect_equal(samples$sample_id, colnames(dset))
    expect_equal(valid_barcodes, rownames(dset))
    expect_equal(as.character(SummarizedExperiment::colData(dset)$minigenes), samples$minigenes)
    expect_true("minigenes" %in% names(SummarizedExperiment::rowData(dset)))

    n_reads = length(readLines(fastq_file)) / 4
    expect_equal(sum(SummarizedExperiment::assay(dset)), n_reads)
    expect_equal(SummarizedExperiment::colData(dset)$total_reads,
                 unname(colSums(SummarizedExperiment::assay(dset))))
    expect_equal(SummarizedExperiment::colData(dset)$mapped_reads,
                 unname(colSums(SummarizedExperiment::assay(dset))))
})

test_that("count_fastq uses read structure orientation", {
    samples = data.frame(sample_id="sample1", minigenes="pat1", rep="1", origin="lib", barcode="GGG")
    constructs = list(pat1 = data.frame(gene_name="gene", mut_id="gene", pep_id="pep",
                                        tiled=TRUE, barcode="TTAACG"))
    fq = tempfile(fileext=".fq")
    writeLines(c("@read1", "GGGNNNNNCGTTAA", "+", "IIIIIIIIIIIIII"), fq)

    parsed = .rs_parse("3B5S6^M")
    expect_false(parsed$sample$revcomp)
    expect_true(parsed$construct$revcomp)

    dset = count_fastq(fq, samples, constructs, "TTAACG", read_structure="3B5S6^M", verbose=FALSE)
    expect_equal(SummarizedExperiment::assay(dset)[1, 1], 1)
})

test_that("count_fastq infers missing read structure", {
    samples = data.frame(sample_id="sample1", minigenes="pat1", rep="1", origin="lib", barcode="GGG")
    constructs = list(pat1 = data.frame(gene_name="gene", mut_id="gene", pep_id="pep",
                                        tiled=TRUE, barcode="TTAACG"))
    fq = tempfile(fileext=".fq")
    writeLines(c("@read1", "GGGNNNNNCGTTAA", "+", "IIIIIIIIIIIIII"), fq)

    dset = count_fastq(fq, samples, constructs, "TTAACG", verbose=FALSE)
    expect_equal(SummarizedExperiment::assay(dset)[1, 1], 1)
})

test_that("count_fastq requires both sample and construct barcode segments", {
    samples = data.frame(sample_id="sample1", minigenes="pat1", rep="1",
                         origin="lib", barcode="GGG")
    constructs = list(pat1=data.frame(gene_name="gene", mut_id="gene", pep_id="pep",
                                     tiled=TRUE, barcode="TTAACG"))
    for (structure in c("3B", "8S6^M"))
        expect_error(count_fastq("unused.fq", samples, constructs, read_structure=structure,
                                 verbose=FALSE), "Counting requires both")
})

test_that("reads with unmatched samples are stored separately by construct or in metadata", {
    samples = data.frame(sample_id="sample1", minigenes="pat1", rep="1",
                         origin="lib", barcode="GGG")
    constructs = list(pat1=data.frame(gene_name=c("gene1", "gene2"),
                                     mut_id=c("gene1", "gene2"), pep_id=c("pep1", "pep2"),
                                     tiled=TRUE, barcode=c("AAAA", "CCCC")))
    fq = tempfile(fileext=".fq")
    on.exit(unlink(fq))
    reads = c("GGGAAAA", "AACAAAA", "TTTAAAA", "AACCCCC", "AACTTTT", "TTTTTTT",
              "GGGTTTT", "GGG")
    writeLines(as.vector(rbind(paste0("@r", seq_along(reads)), reads, "+",
                              strrep("I", nchar(reads)))), fq)

    dset = count_fastq(fq, samples, constructs, valid_barcodes=c("CCCC", "AAAA"),
                       read_structure="3B4M", verbose=FALSE)

    expect_equal(rownames(dset), c("CCCC", "AAAA"))
    expect_equal(SummarizedExperiment::rowData(dset)$unassigned_reads, c(1, 2))
    expect_equal(as.numeric(SummarizedExperiment::assay(dset)), c(0, 1))
    expect_equal(S4Vectors::metadata(dset)$unmapped_reads, 2)
    expect_equal(SummarizedExperiment::colData(dset)$total_reads, 2)
    expect_equal(SummarizedExperiment::colData(dset)$mapped_reads, 1)
})
