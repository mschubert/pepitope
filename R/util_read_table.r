.read_tsv = function(file, ...) {
    utils::read.delim(file, ..., check.names=FALSE, na.strings=c("", "NA"), strip.white=TRUE)
}

.read_samples = function(samples) {
    if (!is.data.frame(samples)) {
        if (is.character(samples) && length(samples) == 1)
            samples = .read_tsv(samples, colClasses=c(sample_id="character", patient="character",
                                                     barcode="character"))
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
        msg = paste("Removing rows" , remove, "because they are empty in the sample sheet")
        warning(msg, immediate.=TRUE)
        samples = samples[!emptyrow,]
    }
    for (name in req)
        if (any(nchar(samples[[name]]) == 0))
            stop(sQuote(name), " must not be empty for any sample in the sample sheet")
    if (any(duplicated(samples$sample_id)))
        stop("'sample_id' values must be unique")
    if (any(duplicated(samples$barcode)))
        stop("'barcode' values must be unique in the sample sheet")

    samples
}
