#' Count barcodes using guide-counter
#'
#' @param tdir  Path to the directory with demultiplexed FASTQ files
#' @param all_constructs  A named list of all construct libraries
#' @param valid_barcodes  A character vector of all possible barcodes
#' @param reverse_complement  Whether to count the reverse complement of the barcodes instead
#' @return      A `SummarizedExperiment` object with counts and metadata
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
count_bc = function(tdir, all_constructs, valid_barcodes, reverse_complement=FALSE) {
    construct_df = merge_constructs(all_constructs)
    if (missing(valid_barcodes))
        valid_barcodes = construct_df$barcode
    if (!is.character(valid_barcodes) && !is.factor(valid_barcodes))
        stop("'valid_barcodes' must be a character vector")
    if (!all(construct_df$barcode %in% valid_barcodes))
        stop("'all_constructs' contains barcodes not in 'valid_barcodes'")

    res = count_external(tdir, valid_barcodes, reverse_complement)

    stats = res$stats |> mutate(sample_id = sub("\\.R1$", "", label))
    counts = data.matrix(res$counts[-(1:2)])
    rownames(counts) = res$counts$guide
    colnames(counts) = sub("\\.R1$", "", colnames(counts))

    samples = readr::read_tsv(file.path(tdir, "samples.tsv"), show_col_types=FALSE)
    meta = samples |>
        left_join(stats |> select(sample_id, total_reads, mapped_reads)) |>
        mutate(patient = factor(patient),
               smp = paste(ifelse(nchar(origin)>8, stringr::word(origin, 1), origin), rep, sep="-"),
               short = paste(patient, smp),
               label = sprintf("%s (%s)", short, sample_id))

    rows = data.frame(barcode=valid_barcodes) |>
        left_join(construct_df) |>
        mutate(bc_type = ifelse(is.na(bc_type), "unused", bc_type),
               bc_type = factor(bc_type, levels=c(names(all_constructs), "unused")))

    SummarizedExperiment(
        list(counts = counts[, meta$sample_id, drop=FALSE]),
        colData = meta,
        rowData = rows
    )
}

#' Bind the tables of all library annotations together
#'
#' @param all_constructs  A named list of all construct libraries
#' @param strict  Whether to make sure the merged table is valid
#' @return  A `data.frame` with all constructs
#'
#' @keywords internal
merge_constructs = function(all_constructs, strict=TRUE) {
    if (!is.list(all_constructs) || is.data.frame(all_constructs))
        stop("'all_constructs' needs to be a named list but is not a list")
    if (is.null(names(all_constructs)))
        stop("'all_constructs' needs to be a named list but has no names")

    expand_barcode = function(name, df) {
        if ("mut_id" %in% colnames(df) && ! "gene_name" %in% colnames(df))
            df$gene_name = sub("_[^_]+$", "", df$mut_id)
        if ("gene_name" %in% colnames(df) && ! "mut_id" %in% colnames(df))
            df$mut_id = df$gene_name
        if ("type" %in% colnames(df) && ! "pep_type" %in% colnames(df))
            df = dplyr::rename(df, pep_type=type)
        req = c("gene_name", "mut_id", "pep_id", "tiled") # pep_type can be NA
        bc_fields = grep("^barcode", colnames(df), value=TRUE)
        if (all(req %in% colnames(df)) && length(bc_fields) > 0) {
            tidyr::pivot_longer(df, bc_fields, values_to="barcode") |> select(-name)
        } else {
            missing = paste(setdiff(req, colnames(df)), collapse=", ")
            if (length(bc_fields) == 0)
                missing = c(missing, "barcode")
            warning("Skipping ", sQuote(name), " because of missing fields: ", missing)
        }
    }
    construct_df = mapply(expand_barcode, name=names(all_constructs),
                          df=all_constructs, SIMPLIFY=FALSE) |>
        bind_rows(.id="bc_type")

    if (strict) {
        if (nrow(construct_df) == 0)
            stop("No data left after merging 'all_constructs'")
        if (any(duplicated(construct_df$barcode)))
            stop("Duplicate barcodes are not allowed but found in 'all_constructs'")
    }

    construct_df
}

#' Use `guide-counter` via a system call to actually count
#'
#' @param tdir  Path to the directory with demultiplexed FASTQ files
#' @param valid_barcodes  A character vector of all possible barcodes
#' @param reverse_complement  Whether to count the reverse complement of the barcodes instead
#' @return      A list with the data.frame meta and matrix counts
#'
#' @importFrom Biostrings reverseComplement DNAStringSet
#' @keywords internal
count_external = function(tdir, valid_barcodes, reverse_complement) {
    tsv = tibble(name = as.character(valid_barcodes)) |>
        mutate(barcode=name, gene=name)
    if (reverse_complement)
        tsv$barcode = as.character(reverseComplement(DNAStringSet(tsv$barcode)))
    lpath = file.path(tdir, "lib.tsv")
    utils::write.table(tsv, file=lpath, sep="\t", row.names=FALSE, quote=FALSE)

    fqs = list.files(tdir, pattern="\\.R1\\.fq\\.gz$", full.names=TRUE)
    if (length(fqs) == 0)
        stop("No fastq files found to count in directory ", sQuote(tdir))

    cmd = paste("guide-counter count",
        "--input", paste(shQuote(fqs), collapse=" "),
        "--library", shQuote(lpath),
        "--offset-min-fraction", shQuote("0.2"),
#        "--exact-match",
        "--output", shQuote(file.path(tdir, "barcodes")))
    system(cmd)

    res = c(counts="barcodes.counts.txt", stats="barcodes.stats.txt") |>
        lapply(\(f) file.path(tdir, f)) |>
        lapply(readr::read_tsv, show_col_types=FALSE)
}
