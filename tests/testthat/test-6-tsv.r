test_that("sample sheets preserve identifiers, metadata names, and quoted text", {
    sheet = tempfile(fileext=".tsv")
    on.exit(unlink(sheet))
    writeLines(c(
        "sample_id\tpatient\trep\torigin\tbarcode\textra field",
        "001\t007\t1\tLibrary\tAAC\t\"O'Brien #1\"",
        "\t\t\t\t\t",
        "002\t008\t2\tMock\tCCC\tNA"
    ), sheet)

    expect_warning(samples <- .read_samples(sheet), "Removing rows 2")
    expect_equal(samples$sample_id, c("001", "002"))
    expect_equal(samples$patient, c("007", "008"))
    expect_equal(samples$rep, c(1, 2))
    expect_equal(samples[["extra field"]], c("O'Brien #1", NA_character_))
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
