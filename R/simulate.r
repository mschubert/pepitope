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
#' @param sample_sheet   The .tsv file containing sample information
#' @param peptide_sheet  The .tsv file containing construct information
#' @param target_reads   How many reads to simulate on average
#' @export
sim_fastq = function(sample_sheet, peptide_sheet, target_reads=1000) {
    sim_reads = function(sample, construct, n) {
        read_seq = paste0(sample$barcode, construct$barcode, construct$tiled)
        Map(paste, sep="\n", USE.NAMES=FALSE,
            paste0("@", sample$barcode, ":", construct$barcode, "_read", seq_len(n)),
            read_seq,
            "+",
            paste0(rep("I", nchar(read_seq)), collapse = "")
        )
    }

    if (missing(sample_sheet))
        sample_sheet = system.file("my_samples.tsv", package="pepitope")
    if (missing(peptide_sheet))
        peptide_sheet = system.file("my_peptides.tsv", package="pepitope")

    samples = readr::read_tsv(sample_sheet, show_col_types=FALSE)
    peptides = readr::read_tsv(peptide_sheet, show_col_types=FALSE)

    idx = expand.grid(p=seq_len(nrow(peptides)), s=seq_len(nrow(samples)))
    idx_s = lapply(idx$s, function(i) samples[i,])
    idx_p = lapply(idx$p, function(i) peptides[i,])
    res = Map(sim_reads, sample=idx_s, construct=idx_p, rpois(nrow(idx), target_reads))

    tdir = tempdir()
    outfile = file.path(tdir, "my_seqdata.fq")
    do.call(cat, c(do.call(c, res), list(file=outfile, sep="\n")))
    outfile
}
