test_that("named sequences stay in their original columns in peptides and reports", {
    subs = GRanges("1", IRanges(1, 1), var_id="1:1_A/T", mut_id="GENE_A2",
                    gene_name="GENE", GENEID="ENSG1", tx_name="tx1",
                    CONSEQUENCE="nonsynonymous", VARAA="V", alt_shift=0)
    subs$ref_nuc = DNAStringSet(c(tx1="ATGGCT"))
    subs$alt_nuc = DNAStringSet(c(tx1="ATGGTT"))
    subs$ref_prot = AAStringSet(c(tx1="MA"))
    subs$alt_prot = AAStringSet(c(tx1="MV"))

    fus = DataFrame(fusion="GENE--OTHER", gene_id_5p="ENSG1", gene_id_3p="ENSG2",
                     tx_id_5p="tx1", tx_id_3p="tx2")
    fus$ref_nuc_5p = DNAStringSet(c(tx1="ATGGCT"))
    fus$ref_nuc_3p = DNAStringSet(c(tx2="ATGTTT"))
    fus$alt_nuc = DNAStringSet(c(fusion="ATGCTT"))

    pep = make_peptides(subs, fus)
    expect_equal(pep$cDNA, c("ATGGCT", "ATGGTT", "ATGGCT", "ATGTTT", "ATGCTT"))
    tiled = pep_tile(pep, tile_size=6, tile_ov=3) |> remove_cutsite(BbsI="GAAGAC")
    report = make_report(subs, subs, fus, tiled)

    expect_false(any(c("ref_nuc", "alt_nuc", "ref_prot", "alt_prot") %in%
                     names(report[["All Variants"]])))
    expect_equal(report[["Unique Protein-Coding"]]$ref_nuc, "ATGGCT")
    expect_equal(report[["Unique Protein-Coding"]]$alt_nuc, "ATGGTT")
    expect_equal(report[["Unique Protein-Coding"]]$ref_prot, "MA")
    expect_equal(report[["Unique Protein-Coding"]]$alt_prot, "MV")
    expect_equal(report[["RNA Fusions"]]$ref_nuc_5p, "ATGGCT")
    expect_equal(report[["RNA Fusions"]]$ref_nuc_3p, "ATGTTT")
    expect_equal(report[["RNA Fusions"]]$alt_nuc, "ATGCTT")
    expect_s4_class(subs$ref_nuc, "DNAStringSet")
    expect_s4_class(fus$ref_nuc_5p, "DNAStringSet")
})
