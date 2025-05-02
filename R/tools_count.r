#' Count barcodes using guide-counter
#'
#' @param tdir  Path to the directory with demultiplexed FASTQ files
#' @param valid_barcodes  A character vector of all possible barcodes
#' @param reverse_complement  Whether to count the reverse complement of the barcodes instead
#' @return      A list with the data.frame meta and matrix counts
count_bc = function(tdir, valid_barcodes, reverse_complement=FALSE) {
    res = count_external(tdir, valid_barcodes, reverse_complement)
    load_dset(res$stats_tsv, res$counts_tsv, res$meta_tsv)
}

#' @keywords internal
count_external = function(tdir, valid_barcodes, reverse_complement) {
    tsv = tibble::tibble(name = as.character(valid_barcodes)) |>
        dplyr::mutate(barcode=name, gene=name)
    if (reverse_complement)
        tsv$barcode = as.character(Biostrings::reverseComplement(DNAStringSet(tsv$barcode)))
    lpath = file.path(tdir, "lib.tsv")
    utils::write.table(tsv, file=lpath, sep="\t", row.names=FALSE, quote=FALSE)

    fqs = list.files(tdir, pattern="\\.fq.gz", full.names=TRUE)

    cmd = paste("guide-counter count",
        "--input", paste(fqs, collapse=" "),
        "--library", lpath,
        "--offset-min-fraction", "0.3",
#        "--exact-match",
        "--output", file.path(tdir, "barcodes"))
    system(cmd)

    res = c(counts="barcodes.counts.txt", stats="barcodes.stats.txt") |>
        lapply(\(f) file.path(tdir, f)) |>
        lapply(readr::read_tsv)
}

#' @keywords internal
load_dset = function(stats_tsv, counts_tsv, meta_tsv) {
    counts_stats = readr::read_tsv(stats_tsv) |> mutate(sample_id = sub("\\.R1$", "", label))
    counts_df = readr::read_tsv(counts_tsv)
    counts = data.matrix(counts_df[-(1:2)])
    rownames(counts) = counts_df$guide
    colnames(counts) = sub("\\.R1$", "", colnames(counts))
    meta = readr::read_tsv(meta_tsv) |>
        left_join(counts_stats |> select(sample_id, total_reads, mapped_reads)) |>
        mutate(patient = factor(patient),
               smp = paste(ifelse(nchar(origin)>8, stringr::word(origin, 1), origin), rep, sep="-"),
               short = paste(patient, smp),
               label = sprintf("%s (%s)", short, sample_id))
    list(meta=meta, counts=counts[,meta$sample_id, drop=FALSE])
}
