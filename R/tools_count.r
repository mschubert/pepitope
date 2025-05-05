#' Count barcodes using guide-counter
#'
#' @param tdir  Path to the directory with demultiplexed FASTQ files
#' @param all_constructs  A named list of all construct libraries
#' @param valid_barcodes  A character vector of all possible barcodes
#' @param reverse_complement  Whether to count the reverse complement of the barcodes instead
#' @return      A list with the data.frame meta and matrix counts
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
count_bc = function(tdir, all_constructs, valid_barcodes, reverse_complement=FALSE) {
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

    construct_df = bind_rows(all_constructs, .id="bc_type")
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

    cmd = paste("guide-counter count",
        "--input", paste(fqs, collapse=" "),
        "--library", lpath,
        "--offset-min-fraction", "0.2",
#        "--exact-match",
        "--output", file.path(tdir, "barcodes"))
    system(cmd)

    res = c(counts="barcodes.counts.txt", stats="barcodes.stats.txt") |>
        lapply(\(f) file.path(tdir, f)) |>
        lapply(readr::read_tsv, show_col_types=FALSE)
}
