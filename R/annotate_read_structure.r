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
#' @export
annotate_read_structure = function(fq, samples, all_constructs, nrec=100000L) {
    samples = .read_samples(samples)
    construct_df = merge_constructs(all_constructs)

    sample_barcodes = toupper(samples$barcode)
    construct_barcodes = toupper(construct_df$barcode)
    .check_barcodes(sample_barcodes, "Sample barcodes")
    .check_barcodes(construct_barcodes, "Construct barcodes")

    reads = readDNAStringSet(fq, format="fastq", nrec=nrec, use.names=FALSE)
    reads = chartr("acgt", "ACGT", reads)
    revcomp = function(seq) as.character(reverseComplement(DNAStringSet(seq)))
    counts = data.frame(
        `B` = .count_positions(reads, sample_barcodes),
        `B<` = .count_positions(reads, revcomp(sample_barcodes)),
        `M` = .count_positions(reads, construct_barcodes),
        `M<` = .count_positions(reads, revcomp(construct_barcodes)),
        check.names=FALSE
    )

    list(
        counts = counts,
        reads = reads,
        structure = .fmt_read_structure(counts, nchar(sample_barcodes[1]), nchar(construct_barcodes[1]))
    )
}

.count_positions = function(reads, barcodes) {
    barcodes = DNAStringSet(barcodes)
    bc_width = BiocGenerics::width(barcodes)[1]
    read_width = BiocGenerics::width(reads)[1]

    sep = DNAString(strrep("N", bc_width - 1L))
    starts = matchPDict(PDict(barcodes), unlist(xscat(reads, sep))) |>
        unlist(use.names=FALSE) |>
        IRanges::start()
    positions = ((starts - 1L) %% (read_width + bc_width - 1L)) + 1L
    counts = tabulate(positions, nbins=read_width)
    counts[(read_width - bc_width + 2L):read_width] = 0L
    counts
}

.best_position = function(counts, op, width) {
    cols = c(op, paste0(op, "<"))
    scores = as.matrix(counts[cols])
    best = arrayInd(which.max(scores), dim(scores))
    if (scores[best] == 0)
        stop("Could not infer read structure because ", op, " barcodes were not found")

    data.frame(start=best[1], width=width, op=cols[best[2]])
}

.fmt_read_structure = function(counts, sample_width, construct_width) {
    sample = .best_position(counts, "B", sample_width)
    construct = .best_position(counts, "M", construct_width)
    segments = rbind(sample, construct) |> dplyr::arrange(start)
    if (any(segments$start[-1] < head(segments$start + segments$width, -1)))
        stop("Could not infer read structure because barcode positions overlap")

    previous_end = c(0L, head(segments$start + segments$width - 1L, -1L))
    skip = segments$start - previous_end - 1L
    skip_tokens = ifelse(skip > 0, paste0(skip, "S"), NA_character_)
    segment_tokens = paste0(segments$width, segments$op)
    tokens = as.vector(rbind(skip_tokens, segment_tokens))
    tokens = stats::na.omit(tokens)
    paste(tokens, collapse="")
}
