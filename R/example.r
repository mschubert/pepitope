#' Create a my_peptides.tsv in the inst directory
#'
#' @keywords internal
example_peptide_file = function() {
    variant_vcf_file = system.file("my_variants.vcf", package="pepitope")
    ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
    asm = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    seqlevelsStyle(ens106) = "UCSC"

    vr1 = readVcfAsVRanges(variant_vcf_file) |>
        filter_variants(min_cov=2, min_af=0.05, pass=TRUE)
    ann = annotate_coding(vr1, ens106, asm)
    subs = ann |> subset_context(15)

    fusion_vcf_file = system.file("my_fusions.vcf", package="pepitope")
    vr2 = readVcfAsVRanges(fusion_vcf_file) |>
        filter_fusions(min_reads=2, min_split_reads=1, min_tools=1)
    seqlevelsStyle(vr2) = "UCSC"
    fus = annotate_fusions(vr2, ens106, asm) |>
        subset_context_fusion(15)

    tiled = make_peptides(subs, fus) |> pep_tile()

    outfile = file.path(system.file(package="pepitope"), "my_peptides.tsv")
    write.table(tiled, outfile, sep="\t", row.names=FALSE, quote=FALSE)
}

#' Create example peptide sheets for multiple samples
#'
#' @param valid_barcodes  A character vector of valid barcodes
#' @return  A named list of peptide/minigene constructs with barcodes
#'
#' @keywords internal
#' @export
example_peptides = function(valid_barcodes, seed=18245) {
    constructs = system.file("my_peptides.tsv", package="pepitope") |>
        readr::read_tsv(show_col_types=FALSE)

    mut_ids = unique(constructs$mut_id)
    if (!is.null(seed))
        set.seed(seed)

    pat1 = constructs |>
        mutate(barcode_1 = valid_barcodes[seq_len(n())],
               barcode_2 = valid_barcodes[seq_len(n()) + n()])

    offset = 2 * nrow(pat1)
    pat2 = constructs |>
        filter(mut_id %in% sample(mut_ids, 15)) |>
        mutate(barcode_1 = valid_barcodes[seq_len(n()) + offset],
               barcode_2 = valid_barcodes[seq_len(n()) + n() + offset])

    offset = offset + 2 * nrow(pat2)
    pat3 = constructs |>
        filter(mut_id %in% sample(mut_ids, 15)) |>
        mutate(barcode_1 = valid_barcodes[seq_len(n()) + offset],
               barcode_2 = valid_barcodes[seq_len(n()) + n() + offset])

    offset = offset + 2 * nrow(pat3)
    common = constructs |>
        filter(mut_id %in% sample(mut_ids, 15)) |>
        mutate(barcode_1 = valid_barcodes[seq_len(n()) + offset],
               barcode_2 = valid_barcodes[seq_len(n()) + n() + offset])

    list(pat1=pat1, pat2=pat2, pat3=pat3, common=common)
}

#' Simulate sequencing data and write them to a temporary FASTQ file
#'
#' @param samples   The .tsv or data.frame file containing sample information
#' @param peptide_sheets  A list, each item containing construct information
#' @param target_reads   How many reads to simulate on average
#' @param custom         Whether to add custom modifications to founds
#' @param seed  The random seed used for sampling the number of reads
#' @return      The path to the created FASTQ file
#'
#' @keywords internal
#' @export
example_fastq = function(samples, peptide_sheets, target_reads=1000, custom=TRUE, seed=91651) {
    sim_reads = function(sample_barcode, construct_seq, n) {
        read_seq = paste0(sample_barcode, construct_seq)
        Map(paste, sep="\n", USE.NAMES=FALSE,
            paste0("@", sample_barcode, ":", construct_seq, "_read", seq_len(n)),
            read_seq,
            "+",
            paste0(rep("I", nchar(read_seq)), collapse = "")
        )
    }

    if (is.character(samples) && length(samples) == 1 && file.exists(samples))
        samples = readr::read_tsv(samples, show_col_types=FALSE)

    all_seq = merge_constructs(peptide_sheets)
    counts = matrix(0, nrow=nrow(all_seq), ncol=length(samples$sample_id),
        dimnames=list(barcode=all_seq$barcode, sample_id=samples$sample_id))

    # fill all expected barcodes (i.e., has an entry in the sample sheet)
    if (!is.null(seed))
        set.seed(seed)
    smp_file = samples |> select(sample_id, patient) |>
        mutate(patient = strsplit(patient, "+", fixed=TRUE)) |>
        tidyr::unnest(patient)
    for (i in seq_len(nrow(smp_file))) {
        matches = all_seq$bc_type == smp_file$patient[i]
        counts[matches, smp_file$sample_id[i]] = rnbinom(sum(matches), mu=target_reads, size=3)
    }

    if (custom) {
        # derive Sample from Mock + noise
        counts[,"screen1"] = counts[,"mock1"] + rnbinom(nrow(counts), mu=counts[,"mock1"], size=5)
        counts[,"screen2"] = counts[,"mock2"] + rnbinom(nrow(counts), mu=counts[,"mock2"], size=5)

        # add specific dropout for NRAS_Q61L in Sample
        bcs = with(all_seq, barcode[bc_type == "pat1" & pep_type == "alt" & mut_id == "NRAS_Q61L"])
        counts[bcs, c("screen1", "screen2")] = rnbinom(4, mu=target_reads/5, size=5)

        # high variance pat2 lib
        pat2_bcs = all_seq$bc_type == "pat2"
        counts[pat2_bcs, "lib1"] = rnbinom(sum(pat2_bcs), mu=target_reads, size=1)

        # lose a quarter of barcodes in pat3 Library and add pat1 contamination
        pat1_bcs = all_seq$bc_type == "pat1"
        pat3_bcs = which(all_seq$bc_type == "pat3") |> head(15)
        counts[pat1_bcs, "lib2"] = rnbinom(sum(pat1_bcs), mu=target_reads/5, size=1)
        counts[pat3_bcs, "lib2"] = rnbinom(length(pat3_bcs), mu=1, size=1)
    }

    # add all non-zero read entries in fastq format
    reads = reshape2::melt(counts) |> filter(value != 0) |>
        inner_join(samples, by=join_by(sample_id))
    res = with(reads, Map(sim_reads, barcode.y, barcode.x, value))

    tdir = tempdir()
    outfile = file.path(tdir, "my_seqdata.fq")
    do.call(cat, c(do.call(c, res), list(file=outfile, sep="\n")))
    outfile
}
