#' Annotate barcode positions in source FASTQ reads
#'
#' @param fq  Path to one FASTQ file
#' @param samples  A sample sheet as `data.frame` in tsv format. Requires the
#'      columns 'sample_id', 'patient', 'rep', 'origin', 'barcode'
#' @param all_constructs  A named list of all construct libraries
#' @return  A `data.frame` with barcode matches per position
#'
#' @importFrom Biostrings readDNAStringSet reverseComplement DNAStringSet subseq
#' @export
annotate_read_structure = function(fq, samples, all_constructs) {
    samples = .read_samples(samples)
    construct_df = merge_constructs(all_constructs)

    sample_barcodes = toupper(samples$barcode)
    construct_barcodes = toupper(construct_df$barcode)
    .check_barcodes(sample_barcodes, "Sample barcodes")
    .check_barcodes(construct_barcodes, "Construct barcodes")

    reads = readDNAStringSet(fq, format="fastq", nrec=100000L, use.names=FALSE)
    reads = chartr("acgt", "ACGT", reads)
    max_width = max(BiocGenerics::width(reads))

    data.frame(
        position = seq_len(max_width),
        sample_fwd = .count_positions(reads, sample_barcodes, max_width),
        sample_rev = .count_positions(reads, .revcomp(sample_barcodes), max_width),
        construct_fwd = .count_positions(reads, construct_barcodes, max_width),
        construct_rev = .count_positions(reads, .revcomp(construct_barcodes), max_width)
    )
}

.count_positions = function(reads, barcodes, max_width) {
    barcodes = DNAStringSet(barcodes)
    width = unique(BiocGenerics::width(barcodes))
    read_widths = BiocGenerics::width(reads)
    positions = seq_len(max_width)
    counts = integer(max_width)
    valid_positions = positions[positions + width - 1L <= max_width]
    counts[valid_positions] = vapply(valid_positions, function(pos) {
        long_enough = read_widths >= pos + width - 1L
        if (!any(long_enough))
            return(0L)
        sum(subseq(reads[long_enough], start=pos, width=width) %in% barcodes)
    }, integer(1))
    counts
}

.revcomp = function(seq) {
    as.character(reverseComplement(DNAStringSet(seq)))
}
