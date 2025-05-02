#' Use fqtk tool to demultiplex fastq files
#'
#' @param fq  A path to the fastq file to demultiplex
#' @param samples  A sample sheet in tsv format
#' @param read_structures  A character string describing the read structure
#' @export
demux_fq = function(fq, samples, read_structures) {
    tdir = tempdir()

    #TODO: create sample sheet if samples is not a file

    cmd = paste("fqtk demux --inputs", fq,
        "--max-mismatches", "0",
        "--read-structures", read_structures,
        "--sample-metadata", samples,
        "--output", tdir)
    system(cmd)

    tdir
}
