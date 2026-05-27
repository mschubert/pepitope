#' Annotate barcode positions in source FASTQ reads
#'
#' @param fq  Path to one FASTQ file
#' @param samples  A sample sheet as `data.frame` in tsv format. Requires the
#'      columns 'sample_id', 'patient', 'rep', 'origin', 'barcode'
#' @param all_constructs  A named list of all construct libraries
#' @param nrec  Number of FASTQ records to inspect
#' @return  A `list` with barcode counts, reads, and inferred read structure
#'
#' @importFrom Biostrings readDNAStringSet reverseComplement DNAStringSet DNAString PDict
#'      matchPDict xscat
#' @importFrom BiocGenerics width
#' @export
annotate_read_structure = function(fq, samples, all_constructs, nrec=100000L) {
    samples = .read_samples(samples)
    construct_df = merge_constructs(all_constructs)
    sample_barcodes = DNAStringSet(toupper(samples$barcode))
    construct_barcodes = DNAStringSet(toupper(construct_df$barcode))
    .check_barcodes(sample_barcodes, "Sample barcodes")
    .check_barcodes(construct_barcodes, "Construct barcodes")

    reads = readDNAStringSet(fq, format="fastq", nrec=nrec, use.names=FALSE)
    reads = chartr("acgt", "ACGT", reads)
    counts = data.frame(
        `B` = .count_positions(reads, sample_barcodes),
        `B<` = .count_positions(reads, reverseComplement(sample_barcodes)),
        `M` = .count_positions(reads, construct_barcodes),
        `M<` = .count_positions(reads, reverseComplement(construct_barcodes)),
        check.names=FALSE
    )

    list(
        counts = counts,
        reads = reads,
        structure = .rs_format(counts, width(sample_barcodes)[1], width(construct_barcodes)[1])
    )
}

.count_positions = function(reads, barcodes) {
    bc_width = width(barcodes)[1]
    read_width = width(reads)[1]

    sep = DNAString(strrep("N", bc_width - 1L))
    starts = matchPDict(PDict(barcodes), unlist(xscat(reads, sep))) |>
        unlist(use.names=FALSE) |>
        IRanges::start()
    positions = ((starts - 1L) %% (read_width + bc_width - 1L)) + 1L
    counts = tabulate(positions, nbins=read_width)
    counts[(read_width - bc_width + 2L):read_width] = 0L
    counts
}
