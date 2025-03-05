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
