# provide fastq file
# use fqtk cli to demux fastq files
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
