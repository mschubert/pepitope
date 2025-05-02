#' Create a peptide .tsv in the inst directory
#'
#' @param outfile
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

    bcs = make_lib() |> pull(barcode) |> head(nrow(tiled))
    res = tiled |> dplyr::select(-cDNA) |> dplyr::mutate(barcode = bcs)
    write.table(res, outfile, sep="\t", row.names=FALSE, quote=FALSE)
}

#' Simulate sequencing data and write them to a FASTQ file
#'
#' @param outfile        The .fastq.gz file to save results to
#' @param sample_sheet   The .tsv file containing sample information
#' @param peptide_sheet  The .tsv file containing construct information
#' @param target_reads   How many reads to simulate on average
#' @keywords internal
sim_fastq = function(outfile, sample_sheet, peptide_sheet, target_reads=100) {
    if (missing(sample_sheet))
        sample_sheet = system.file("my_samples.tsv", package="pepitope")
    if (missing(peptide_sheet))
        peptide_sheet = system.file("my_peptides.tsv", package="pepitope")

    samples = readr::read_tsv(sample_sheet)
    peptides = readr::read_tsv(peptide_sheet)

    # multiple samples: noise + 2 clear drop-outs

    # merge with barcodes

    # use biostrings to write fastq
}
