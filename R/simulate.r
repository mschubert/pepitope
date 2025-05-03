#' Create a peptide .tsv in the inst directory
#'
#' @param outfile  The file to save peptide table to
#' @keywords internal
make_pep_table = function(outfile) {
    if (missing(outfile))
        outfile = file.path(system.file(package="pepitope"), "my_peptides.tsv")

    variant_vcf_file = system.file("my_variants.vcf", package="pepitope")
    ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
    asm = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    seqlevelsStyle(ens106) = "UCSC"

    vr1 = readVcfAsVRanges(variant_vcf_file) |>
        filter_variants(min_cov=2, min_af=0.05, pass=TRUE)
    ann = annotate_coding(vr1, ens106, asm)
    subs = ann |> subset_context(15)
    tiled = pep_tile(subs) |> remove_cutsite(BbsI="GAAGAC")

    bc_file = "https://raw.githubusercontent.com/hawkjo/freebarcodes/master/barcodes/barcodes12-1.txt"
    bcs = readr::read_tsv(bc_file, col_names=FALSE, show_col_types=FALSE)$X1 |> head(nrow(tiled))
    res = tiled |> dplyr::select(-cDNA) |> dplyr::mutate(barcode = bcs)
    write.table(res, outfile, sep="\t", row.names=FALSE, quote=FALSE)
}

#' Simulate sequencing data and write them to a FASTQ file
#'
#' @param sample_sheet   The .tsv or data.frame file containing sample information
#' @param peptide_sheet  A list, each item containing construct information
#' @param target_reads   How many reads to simulate on average
#' @export
sim_fastq = function(samples, peptide_sheets, target_reads=1000) {
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

    all_seq = dplyr::bind_rows(peptide_sheets, .id="bc_type")
    counts = matrix(0, nrow=nrow(all_seq), ncol=length(samples$sample_id),
        dimnames=list(barcode=all_seq$barcode, sample_id=samples$sample_id))

    # fill all expected barcodes (i.e., has an entry in the sample sheet)
    smp_file = samples |> select(sample_id, patient) |>
        mutate(patient = strsplit(patient, "+", fixed=TRUE)) |>
        tidyr::unnest(patient)
    for (i in seq_len(nrow(smp_file))) {
        matches = all_seq$bc_type == smp_file$patient[i]
        counts[matches, smp_file$sample_id[i]] = rpois(sum(matches), target_reads)
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
