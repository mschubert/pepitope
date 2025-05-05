#' Use fqtk tool to demultiplex fastq files
#'
#' @param fq  A path to the fastq file to demultiplex
#' @param samples  A sample sheet as `data.frame` in tsv format
#' @param read_structures  A character string describing the read structure
#' @export
demux_fq = function(fq, samples, read_structures) {
    tdir = tempdir()

    fname = file.path(tdir, "samples.tsv")
    if (is.character(samples) && length(samples) == 1) {
        file.copy(samples, fname)
    } else if (is.data.frame(samples)) {
        utils::write.table(samples, file=fname, row.names=FALSE, sep="\t")
        samples = fname
    } else {
        stop("'samples' argument needs to be a single tsv file or a data.frame")
    }

    cmd = paste("fqtk demux --inputs", fq,
        "--max-mismatches", "0",
        "--read-structures", read_structures,
        "--sample-metadata", samples,
        "--output", tdir)
    system(cmd)

    tdir
}
