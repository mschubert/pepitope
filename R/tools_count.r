# start with demuxed fastq files in temp dir
# use guide_counter to count barcodes
count_bc = function(tdir, lib) {
    lpath = file.path(tdir, "lib.tsv")
    utils::write.table(lib, file=lpath, sep="\t", row.names=FALSE, quote=FALSE)
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

load_dset = function(args) {
    counts_stats = readr::read_tsv(args$stats_tsv) |> mutate(sample_id = sub("\\.R1$", "", label))
    counts_df = readr::read_tsv(args$counts_tsv)
    counts = data.matrix(counts_df[-(1:2)])
    rownames(counts) = counts_df$guide
    colnames(counts) = sub("\\.R1$", "", colnames(counts))
    meta = readr::read_tsv(args$meta_tsv) |>
        left_join(counts_stats |> select(sample_id, total_reads, mapped_reads)) |>
        mutate(patient = factor(patient),
               smp = paste(ifelse(nchar(origin)>8, stringr::word(origin, 1), origin), rep, sep="-"),
               short = paste(patient, smp),
               label = sprintf("%s (%s)", short, sample_id))
    list(meta=meta, counts=counts[,meta$sample_id, drop=FALSE])
}
