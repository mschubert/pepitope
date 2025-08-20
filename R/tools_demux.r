#' Use fqtk tool to demultiplex fastq files
#'
#' @param fq  A path to the fastq file to demultiplex
#' @param samples  A sample sheet as `data.frame` in tsv format. Requires the
#'      columns 'sample_id', 'patient', 'rep', 'origin', 'barcode'
#' @param read_structures  A character string describing the read structure
#'
#' @export
demux_fq = function(fq, samples, read_structures) {
    tdir = tempfile()
    dir.create(tdir)

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
        warning(paste("Removing rows" , remove, "because they are empty in the sample sheet"))
        samples = samples[!emptyrow,]
    }
    for (name in req)
        if (any(nchar(samples[[name]]) == 0))
            stop(sQuote(name), " must not be empty for any sample in the sample sheet")

    sample_tsv = file.path(tdir, "samples.tsv")
    utils::write.table(samples, file=sample_tsv, row.names=FALSE, sep="\t")

    fqtkWrapper::fqtk_demux(
        inputs = fq,
        max_mismatches = 0L,
        read_structures = read_structures,
        sample_metadata = sample_tsv,
        output = tdir
    )
    tdir
}
