#' Annotate barcode positions in source FASTQ reads
#'
#' If structure inference fails, a warning is emitted and the structure is `"unknown"`.
#'
#' @param fq  Path to one FASTQ file
#' @param sample_barcodes  A character vector of known sample barcodes
#' @param construct_barcodes  A character vector of known construct barcodes
#' @param nrec  Number of FASTQ records to inspect
#' @param min_frac  Minimum fraction of inspected reads matching each barcode type
#'      at its inferred position in either orientation
#' @return  A `list` with barcode counts, reads, and inferred read structure
#'
#' @importFrom Biostrings readDNAStringSet reverseComplement DNAStringSet DNAString
#' @importFrom Biostrings PDict matchPDict xscat
#' @importFrom BiocGenerics width
#' @keywords internal
.rs_annotate = function(fq, sample_barcodes, construct_barcodes, nrec=100000L, min_frac=0.2) {
    sample_barcodes = DNAStringSet(toupper(sample_barcodes))
    construct_barcodes = DNAStringSet(toupper(construct_barcodes))
    .check_barcodes(sample_barcodes, "Sample barcodes")
    .check_barcodes(construct_barcodes, "Construct barcodes")

    reads = readDNAStringSet(fq, format="fastq", nrec=nrec, use.names=FALSE)
    reads = chartr("acgt", "ACGT", reads)
    counts = data.frame(
        `B` = .rs_count(reads, sample_barcodes),
        `^B` = .rs_count(reads, reverseComplement(sample_barcodes)),
        `M` = .rs_count(reads, construct_barcodes),
        `^M` = .rs_count(reads, reverseComplement(construct_barcodes)),
        check.names=FALSE
    )

    structure = tryCatch(
        .rs_format(counts, width(sample_barcodes)[1], width(construct_barcodes)[1],
                   n_reads=length(reads), min_frac=min_frac),
        error=function(e) { warning(conditionMessage(e), immediate.=TRUE); "unknown" }
    )
    list(counts=counts, reads=reads, structure=structure)
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

#' Match inferred barcode segments in individual reads
#'
#' @param reads  Read sequences as a character vector or `DNAStringSet`
#' @param structure  An inferred read structure string, including `"unknown"`
#' @param sample_barcodes  Known sample barcodes
#' @param construct_barcodes  Known construct barcodes
#' @return  A `data.frame` with read indices, barcode types, ranges, orientations, and
#'      logical `present` values; unknown structures return zero rows
#'
#' @importFrom Biostrings DNAStringSet subseq reverseComplement
#' @keywords internal
.rs_match = function(reads, structure, sample_barcodes, construct_barcodes) {
    if (structure == "unknown")
        return(data.frame(read=integer(), start=integer(), width=integer(),
                          revcomp=logical(), type=character(), present=logical()))

    structure = .rs_parse(structure)
    features = rbind(data.frame(structure$sample, type=rep("B", nrow(structure$sample))),
                     data.frame(structure$construct, type=rep("M", nrow(structure$construct))))
    matches = merge(data.frame(read=seq_along(reads)), features)
    reads = DNAStringSet(toupper(as.character(reads)))
    barcodes = subseq(reads[matches$read], start=matches$start, width=matches$width)
    barcodes[matches$revcomp] = reverseComplement(barcodes[matches$revcomp])
    matches$present = ifelse(matches$type == "B",
                             as.character(barcodes) %in% toupper(sample_barcodes),
                             as.character(barcodes) %in% toupper(construct_barcodes))
    matches
}

#' Parse a read structure
#'
#' @param read_structures  A character string describing the FASTQ read structure
#' @return  A `list` with sample and construct barcode ranges; missing types have zero rows
#'
#' @keywords internal
.rs_parse = function(read_structures) {
    if (!is.character(read_structures) || length(read_structures) != 1 || is.na(read_structures))
        stop("'read_structures' must be a single character string")

    operators = c("T", "B", "M", "C", "S")
    matches = gregexpr("([0-9]+|\\+)(\\^)?([TBMCS])", read_structures, perl=TRUE)[[1]]
    if (matches[1] == -1)
        stop("Invalid read structure: ", sQuote(read_structures))
    tokens = regmatches(read_structures, list(matches))[[1]]
    if (paste(tokens, collapse="") != read_structures)
        stop("Invalid read structure: ", sQuote(read_structures))

    starts = attr(matches, "capture.start")
    lengths = attr(matches, "capture.length")
    ranges = lapply(seq_along(tokens), function(i) {
        len_text = substr(read_structures, starts[i, 1], starts[i, 1] + lengths[i, 1] - 1L)
        revcomp = lengths[i, 2] > 0
        op = substr(read_structures, starts[i, 3], starts[i, 3] + lengths[i, 3] - 1L)
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
        data.frame(width=len, op=op, revcomp=revcomp)
    })

    ranges = do.call(rbind, ranges)
    ranges$start = c(1L, utils::head(cumsum(ranges$width) + 1L, -1L))
    ranges = ranges[c("start", "width", "op", "revcomp")]
    sample = ranges[ranges$op == "B", c("start", "width", "revcomp")]
    construct = ranges[ranges$op == "M", c("start", "width", "revcomp")]
    if (nrow(sample) > 1)
        stop("'read_structures' must contain at most one sample barcode ('B') segment")
    if (nrow(construct) > 1)
        stop("'read_structures' must contain at most one construct barcode ('M') segment")
    if (nrow(sample) + nrow(construct) == 0)
        stop("'read_structures' must contain a sample ('B') or construct ('M') barcode segment")
    list(sample=sample, construct=construct)
}

#' Format a read structure from barcode position counts
#'
#' Raises an error if neither barcode type meets `min_frac` or inferred barcode positions overlap.
#'
#' @param counts  A `data.frame` with `B`, `^B`, `M`, and `^M` count columns
#' @param sample_width  Width of the sample barcode
#' @param construct_width  Width of the construct barcode
#' @param n_reads  Number of inspected reads
#' @param min_frac  Minimum fraction of reads matching each barcode type at its
#'      inferred position in either orientation
#' @return  A read structure string containing the barcode types that meet `min_frac`
#'
#' @keywords internal
.rs_format = function(counts, sample_width, construct_width, n_reads, min_frac=0.2) {
    best_position = function(op, width) {
        cols = c(op, paste0("^", op))
        scores = as.matrix(counts[cols])
        best = arrayInd(which.max(scores), dim(scores))
        data.frame(start=best[1], width=width, op=cols[best[2]], frac=scores[best] / n_reads)
    }

    sample = best_position("B", sample_width)
    construct = best_position("M", construct_width)
    segments = rbind(sample, construct)
    segments = segments[segments$frac >= min_frac, ]
    if (nrow(segments) == 0)
        stop("Could not infer read structure: barcode match fraction below min_frac (", min_frac, ")")
    segments = dplyr::arrange(segments, start)
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
