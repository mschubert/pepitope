#TODO:
# * write system calls for these tools
#   * include rust crates with extendr
#   * call the api instead of the cli

get_phased_variants = function() {
    # provide .bam, .gtf files
    # use microphaser cli to get phased variants
}

fq = "/data/groups/wg-schubert/shared/Endo+Penile_data/7977_1_Oligolibrary_QC_endo_penis_TAA_TGCCGTTACT-GCCATATATC_S18_R1_001.fastq.gz"
samples = "2024-12_plasmid-Endo.tsv"
readr::read_tsv("2024-12_plasmid-Endo.tsv")

# make barcode lib
# https://github.com/hawkjo/freebarcodes
# https://www.pnas.org/doi/full/10.1073/pnas.1802640115
make_lib = function() {
    url = "https://raw.githubusercontent.com/hawkjo/freebarcodes/master/barcodes/barcodes12-1.txt"
    bcs = read.table(url)$V1

    tsv = tibble::tibble(oligo_id = paste0("oligo_", seq_along(bcs)),
                 barcode = as.character(bcs)) |>
        dplyr::mutate(gene = oligo_id)
}

# provide fastq file
# use fqtk cli to demux fastq files
demux_fq = function(fq, samples) {
    tdir = tempdir()

    #TODO: create sample sheet if samples is not a file

    cmd = paste("fqtk demux --inputs", fq,
        "--max-mismatches", "0",
        "--read-structures", "7B+T",
        "--sample-metadata", samples,
        "--output", tdir)
    system(cmd)

    tdir
}

# start with demuxed fastq files in temp dir
# use guide_counter to count barcodes
count_bc = function(tdir, lib) {
    lpath = file.path(tdir, "lib.tsv")
    write.table(lib, file=lpath, sep="\t", row.names=FALSE, quote=FALSE)
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
