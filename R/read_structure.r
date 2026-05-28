#' Annotate barcode positions in source FASTQ reads
#'
#' @param fq  Path to one FASTQ file
#' @param samples  A sample sheet as `data.frame` in tsv format. Requires the
#'      columns 'sample_id', 'patient', 'rep', 'origin', 'barcode'
#' @param all_constructs  A named list of all construct libraries
#' @param nrec  Number of FASTQ records to inspect
#' @return  A `list` with barcode counts, reads, and inferred read structure
#'
#' @importFrom Biostrings readDNAStringSet reverseComplement DNAStringSet DNAString PDict matchPDict xscat
#' @importFrom BiocGenerics width
#' @keywords internal
.rs_annotate = function(fq, samples, all_constructs, nrec=100000L) {
    samples = .read_samples(samples)
    construct_df = merge_constructs(all_constructs)
    sample_barcodes = DNAStringSet(toupper(samples$barcode))
    construct_barcodes = DNAStringSet(toupper(construct_df$barcode))
    .check_barcodes(sample_barcodes, "Sample barcodes")
    .check_barcodes(construct_barcodes, "Construct barcodes")

    reads = readDNAStringSet(fq, format="fastq", nrec=nrec, use.names=FALSE)
    reads = chartr("acgt", "ACGT", reads)
    counts = data.frame(
        `B` = .rs_count(reads, sample_barcodes),
        `B<` = .rs_count(reads, reverseComplement(sample_barcodes)),
        `M` = .rs_count(reads, construct_barcodes),
        `M<` = .rs_count(reads, reverseComplement(construct_barcodes)),
        check.names=FALSE
    )

    list(
        counts = counts,
        reads = reads,
        structure = .rs_format(counts, width(sample_barcodes)[1], width(construct_barcodes)[1])
    )
}

.rs_count = function(reads, barcodes) {
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

#' Parse a read structure
#'
#' @param read_structures  A character string describing the FASTQ read structure
#' @return  A `list` with sample and construct barcode ranges
#'
#' @keywords internal
.rs_parse = function(read_structures) {
    if (!is.character(read_structures) || length(read_structures) != 1 || is.na(read_structures))
        stop("'read_structures' must be a single character string")

    operators = c("T", "B", "M", "C", "S")
    matches = gregexpr("([0-9]+|\\+)([TBMCS])(<)?", read_structures, perl=TRUE)[[1]]
    if (matches[1] == -1)
        stop("Invalid read structure: ", sQuote(read_structures))
    tokens = regmatches(read_structures, list(matches))[[1]]
    if (paste(tokens, collapse="") != read_structures)
        stop("Invalid read structure: ", sQuote(read_structures))

    starts = attr(matches, "capture.start")
    lengths = attr(matches, "capture.length")
    pos = 1L
    ranges = lapply(seq_along(tokens), function(i) {
        len_text = substr(read_structures, starts[i, 1], starts[i, 1] + lengths[i, 1] - 1L)
        op = substr(read_structures, starts[i, 2], starts[i, 2] + lengths[i, 2] - 1L)
        revcomp = lengths[i, 3] > 0
        if (!op %in% operators)
            stop("Unsupported read structure operator: ", sQuote(op))
        if (len_text == "+") {
            if (i != length(tokens))
                stop("'+' may only be used in the final read structure token")
            if (op %in% c("B", "M"))
                stop("'+' cannot be used for barcode operators 'B' or 'M'")
            len = NA_integer_
        } else {
            len = as.integer(len_text)
            if (is.na(len) || len < 1)
                stop("Read structure lengths must be positive integers")
        }
        res = data.frame(start=pos, width=len, op=op, revcomp=revcomp)
        pos <<- if (is.na(len)) NA_integer_ else pos + len
        res
    })

    ranges = do.call(rbind, ranges)
    sample = ranges[ranges$op == "B", c("start", "width", "revcomp")]
    construct = ranges[ranges$op == "M", c("start", "width", "revcomp")]
    if (nrow(sample) != 1)
        stop("'read_structures' must contain exactly one sample barcode ('B') segment")
    if (nrow(construct) != 1)
        stop("'read_structures' must contain exactly one construct barcode ('M') segment")
    list(sample=sample, construct=construct)
}

#' Format a read structure from barcode position counts
#'
#' @param counts  A `data.frame` with `B`, `B<`, `M`, and `M<` count columns
#' @param sample_width  Width of the sample barcode
#' @param construct_width  Width of the construct barcode
#' @return  A read structure string
#'
#' @keywords internal
.rs_format = function(counts, sample_width, construct_width) {
    best_position = function(op, width) {
        cols = c(op, paste0(op, "<"))
        scores = as.matrix(counts[cols])
        best = arrayInd(which.max(scores), dim(scores))
        if (scores[best] == 0)
            stop("Could not infer read structure because ", op, " barcodes were not found")

        data.frame(start=best[1], width=width, op=cols[best[2]])
    }

    sample = best_position("B", sample_width)
    construct = best_position("M", construct_width)
    segments = rbind(sample, construct) |> dplyr::arrange(start)
    if (any(segments$start[-1] < utils::head(segments$start + segments$width, -1)))
        stop("Could not infer read structure because barcode positions overlap")

    previous_end = c(0L, utils::head(segments$start + segments$width - 1L, -1L))
    skip = segments$start - previous_end - 1L
    skip_tokens = ifelse(skip > 0, paste0(skip, "S"), NA_character_)
    segment_tokens = paste0(segments$width, segments$op)
    tokens = as.vector(rbind(skip_tokens, segment_tokens))
    tokens = stats::na.omit(tokens)
    paste(tokens, collapse="")
}
