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

    tiled = pep_tile(pep, tile_size=6, tile_ov=3) |> remove_cutsite(EcoRI="GAATTC")
    report = make_report(subs, subs, fus, tiled)
    expect_true("EcoRI_replaced" %in% names(report[["93 nt Peptides"]]))
    expect_false("BbsI_replaced" %in% names(report[["93 nt Peptides"]]))
})

test_that("cut site removal achieves the minimum number of nucleotide substitutions", {
    # Exhaustive enumeration is an independent reference for these small examples.
    code = getGeneticCode()
    cases = list(
        c("GAAGAC", "GAAGAC"),
        c("GTCTTC", "GAAGAC"),
        c("GAAGACGTCTTC", "GAAGAC"),
        c("AGAAGACGAAGACAAA", "GAAGAC"),
        c("AGAAGACGACAA", "GAAGAC"),
        c("AAGTCTTCTTCA", "GAAGAC"),
        c("AAAAAA", "AAAA"),
        c("AAAAAAA", "AAAA"),
        c("GAAGACGAATTC", "GAAGAC", "GAATTC"),
        c("GAAGAAGAAGAA", "GAAGAA", "GAGGAA"),
        c("AAAAAA", "AAAA", "AAA")
    )
    for (case in cases) {
        nuc = case[1]
        site = case[-1]
        starts = seq_len(nchar(nuc) %/% 3L) * 3L - 2L
        codons = substring(nuc, starts, starts + 2L)
        choices = lapply(codons, function(codon) names(code)[code == code[[codon]]])
        candidates = do.call(paste0, do.call(expand.grid, c(choices, stringsAsFactors=FALSE)))
        candidates = paste0(candidates, substring(nuc, length(codons) * 3L + 1L))
        patterns = unique(c(site, as.character(reverseComplement(DNAStringSet(site)))))
        valid = Reduce(`&`, lapply(patterns, function(s) !grepl(s, candidates, fixed=TRUE)))
        bases = strsplit(nuc, "")[[1]]
        distance = function(seq) sum(strsplit(seq, "")[[1]] != bases)
        optimum = min(vapply(candidates[valid], distance, integer(1)))

        for (seed in c(0, 1, 178529)) {
            result = remove_cutsite_nuc(nuc, site, seed=seed)
            expect_true(result %in% candidates[valid], info=nuc)
            expect_equal(distance(result), optimum, info=nuc)
        }
    }
})

test_that("cut site removal avoids newly created sites and overlapping occurrences", {
    # A synonymous CGA -> AGA edit moves the BbsI site from position 2 to position 5.
    original = "AGAAGACGACAA"
    moved = "AGAAGAAGACAA"
    expect_equal(translate(DNAString(original)), translate(DNAString(moved)))
    expect_equal(IRanges::start(Biostrings::matchPattern("GAAGAC", DNAString(original))), 2L)
    expect_equal(IRanges::start(Biostrings::matchPattern("GAAGAC", DNAString(moved))), 5L)
    result = remove_cutsite_nuc(original, "GAAGAC")
    expect_equal(vcountPattern("GAAGAC", result) + vcountPattern("GTCTTC", result), 0L)

    # Three overlapping AAAA sites can all be broken by a single synonymous change.
    result = remove_cutsite_nuc("AAAAAA", "AAAA")
    expect_equal(vcountPattern("AAAA", result), 0L)
    expect_equal(sum(strsplit(result, "")[[1]] != rep("A", 6)), 1L)
})

test_that("cut site removal batches sequences and preserves original site counts", {
    sequences = c("ATGGCCGCC", "GAAGAC", "GTCTTC", "GAAGACGTCTTC")
    pep = data.frame(pep_id=paste0("p", seq_along(sequences)), tiled=sequences)
    result = remove_cutsite(pep, BbsI="GAAGAC")
    expect_equal(result$BbsI_replaced, c(0L, 1L, 1L, 2L))
    expect_equal(result$tiled[1], sequences[1])
    expect_equal(translate(DNAStringSet(result$tiled), no.init.codon=TRUE),
                 translate(DNAStringSet(sequences), no.init.codon=TRUE))
    expect_equal(result$tiled,
                 vapply(sequences, remove_cutsite_nuc, character(1), site="GAAGAC",
                        seed=178529, USE.NAMES=FALSE))
    expect_equal(vcountPattern("GAAGAC", result$tiled), rep(0L, 4))
    expect_equal(vcountPattern("GTCTTC", result$tiled), rep(0L, 4))
    expect_equal(remove_cutsite(pep[FALSE, ], BbsI="GAAGAC")$tiled, character())
})

test_that("many cut sites are removed reproducibly without altering the R random state", {
    sequence = strrep("GAAGACGTCTTC", 1000)
    set.seed(123)
    random_state = .Random.seed
    result = remove_cutsite_nuc(sequence, "GAAGAC", seed=178529)
    expect_identical(.Random.seed, random_state)
    expect_identical(remove_cutsite_nuc(sequence, "GAAGAC", seed=178529), result)
    expect_equal(vcountPattern("GAAGAC", result) + vcountPattern("GTCTTC", result), 0L)
    expect_equal(translate(DNAString(result)), translate(DNAString(sequence)))
    expect_equal(sum(strsplit(result, "")[[1]] != strsplit(sequence, "")[[1]]), 2000L)
})

test_that("named cut sites are removed together and counted separately on both strands", {
    sequences = c("GAAGACGAATTCGGTCTC", "GTCTTCGAGACC", "ATGGCCGCC")
    pep = data.frame(pep_id=paste0("p", seq_along(sequences)), tiled=sequences)
    sites = c(BbsI="GAAGAC", EcoRI="GAATTC", BsaI="GGTCTC")
    result = remove_cutsite(pep, BbsI=sites[[1]], EcoRI=sites[[2]], BsaI=sites[[3]])
    expect_equal(result$BbsI_replaced, c(1L, 1L, 0L))
    expect_equal(result$EcoRI_replaced, c(1L, 0L, 0L))
    expect_equal(result$BsaI_replaced, c(1L, 1L, 0L))
    expect_equal(translate(DNAStringSet(result$tiled), no.init.codon=TRUE),
                 translate(DNAStringSet(sequences), no.init.codon=TRUE))
    for (site in unique(c(sites, as.character(reverseComplement(DNAStringSet(sites))))))
        expect_equal(vcountPattern(site, result$tiled), rep(0L, 3))

    # Removing GAAGAA must not create GAGGAA, which is absent from the original.
    pep = data.frame(pep_id="p1", tiled="GAAGAAGAAGAA")
    result = remove_cutsite(pep, Repeat="GAAGAA", Created="GAGGAA")
    expect_equal(result$Repeat_replaced, 3L)
    expect_equal(result$Created_replaced, 0L)
    expect_equal(sum(strsplit(result$tiled, "")[[1]] != strsplit(pep$tiled, "")[[1]]), 3L)
})

test_that("cut site removal reports impossible designs and preserves incomplete codons", {
    expect_error(remove_cutsite_nuc("ATG", "ATG"), "No synonymous sequence")
    expect_error(remove_cutsite_nuc("ATGAA", "AA"), "No synonymous sequence")
    result = remove_cutsite_nuc("AAAAAAA", "AAAA")
    expect_equal(substring(result, 7), "A")
    expect_equal(nchar(result), 7L)
    expect_error(remove_cutsite_nuc("AAAAAA", ""), "recognition sequence")
    pep = data.frame(pep_id="p1", tiled="GAAGAC")
    expect_error(remove_cutsite(pep), "unique, nonempty names")
    expect_error(remove_cutsite(pep, "GAAGAC"), "unique, nonempty names")
    expect_error(remove_cutsite(pep, Site="GAAGAC", Site="GAATTC"), "unique, nonempty names")
})
