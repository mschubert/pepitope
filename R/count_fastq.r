#' Count barcodes directly from source FASTQ files
#'
#' @param fq  Path to one FASTQ file
#' @param samples  A sample sheet as `data.frame` in tsv format. Requires the
#'      columns 'sample_id', 'patient', 'rep', 'origin', 'barcode'
#' @param all_constructs  A named list of all construct libraries
#' @param valid_barcodes  A character vector of all possible construct barcodes
#' @param read_structures  A character string describing the FASTQ read structure. `B` segments
#'      are matched against sample barcodes and `M` segments against construct barcodes.
#' @param reverse_complement  Whether to count the reverse complement of the construct barcodes instead
#' @param verbose  Whether to print progress messages (default: TRUE)
#' @return      A `SummarizedExperiment` object with counts and metadata
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Biostrings reverseComplement DNAStringSet
#' @export
count_fastq = function(fq, samples, all_constructs, valid_barcodes, read_structures="7B12M+T",
                       reverse_complement=FALSE, verbose=TRUE) {
    samples = .read_samples(samples)
    all_samples = strsplit(samples$patient, "+", fixed=TRUE) |> unlist()
    missing = setdiff(all_samples, names(all_constructs))
    if (length(missing) > 0)
        stop("Missing minigene annotations for: ", paste(sQuote(missing), collapse=", "))

    construct_df = merge_constructs(all_constructs)
    if (missing(valid_barcodes))
        valid_barcodes = construct_df$barcode
    if (!is.character(valid_barcodes) && !is.factor(valid_barcodes))
        stop("'valid_barcodes' must be a character vector")
    if (!all(construct_df$barcode %in% valid_barcodes))
        stop("'all_constructs' contains barcodes not in 'valid_barcodes'")
    if (!is.character(fq) || length(fq) != 1)
        stop("'fq' argument needs to be a single FASTQ file")
    fq = path.expand(fq)

    read_structure = .rs_parse(read_structures)

    count_barcodes = as.character(valid_barcodes)
    if (reverse_complement)
        count_barcodes = as.character(reverseComplement(DNAStringSet(count_barcodes)))
    sample_barcodes = toupper(samples$barcode)
    count_barcodes = toupper(count_barcodes)
    .check_barcodes(sample_barcodes, "Sample barcodes")
    .check_barcodes(count_barcodes, "Construct barcodes")
    if (any(nchar(sample_barcodes) != sum(read_structure$sample$width)))
        stop("Sample barcode width does not match 'B' segments in 'read_structures'")
    if (any(nchar(count_barcodes) != sum(read_structure$construct$width)))
        stop("Construct barcode width does not match 'M' segments in 'read_structures'")

    res = count_fastq_barcodes_cpp(
        fq = fq,
        sample_barcodes = stats::setNames(sample_barcodes, samples$sample_id),
        construct_barcodes = stats::setNames(count_barcodes, as.character(valid_barcodes)),
        sample_start = read_structure$sample$start,
        sample_width = read_structure$sample$width,
        construct_start = read_structure$construct$start,
        construct_width = read_structure$construct$width,
        verbose = verbose
    )

    counts = res$counts
    rownames(counts) = valid_barcodes
    colnames(counts) = samples$sample_id

    stats = data.frame(sample_id = samples$sample_id,
                       total_reads = res$total_reads,
                       mapped_reads = res$mapped_reads)
    meta = samples |>
        left_join(stats, by="sample_id") |>
        mutate(patient = factor(patient),
               smp = paste(ifelse(nchar(origin)>8, stringr::word(origin, 1), origin), rep, sep="-"),
               short = paste(patient, smp),
               label = sprintf("%s (%s)", short, sample_id))

    rows = data.frame(barcode=valid_barcodes) |>
        left_join(construct_df, by="barcode") |>
        mutate(bc_type = ifelse(is.na(bc_type), "unused", bc_type),
               bc_type = factor(bc_type, levels=c(names(all_constructs), "unused")))

    SummarizedExperiment(
        list(counts = counts[, meta$sample_id, drop=FALSE]),
        colData = meta,
        rowData = rows
    )
}

.read_samples = function(samples) {
    if (!is.data.frame(samples)) {
        if (is.character(samples) && length(samples) == 1)
            samples = readr::read_tsv(samples, show_col_types=FALSE)
        else
            stop("'samples' argument needs to be a single tsv file or a data.frame")
    }

    req = c("sample_id", "patient", "rep", "origin", "barcode")
    missing = setdiff(req, colnames(samples))
    if (length(missing) > 0)
        stop("Required columns not found in sample sheet: ", paste(sQuote(missing), collapse=", "))
    emptyrow = apply(samples, 1, function(r) all(is.na(r)))
    if (any(emptyrow)) {
        remove = paste(which(emptyrow), collapse=", ")
        msg = paste("Removing rows" , remove, "because they are empty in the sample sheet")
        warning(msg, immediate.=TRUE)
        samples = samples[!emptyrow,]
    }
    for (name in req)
        if (any(nchar(samples[[name]]) == 0))
            stop(sQuote(name), " must not be empty for any sample in the sample sheet")
    if (any(duplicated(samples$sample_id)))
        stop("'sample_id' values must be unique")
    if (any(duplicated(samples$barcode)))
        stop("'barcode' values must be unique in the sample sheet")

    samples
}

.check_barcodes = function(barcodes, label) {
    if (length(barcodes) == 0)
        stop(label, " must not be empty")
    if (anyNA(barcodes))
        stop(label, " must not contain NA values")
    widths = nchar(barcodes)
    if (any(widths == 0))
        stop(label, " must not contain empty strings")
    if (length(unique(widths)) != 1)
        stop(label, " must all have the same width")
    if (any(duplicated(barcodes)))
        stop(label, " must not contain duplicate values")
}
