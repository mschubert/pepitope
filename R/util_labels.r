.format_sample_labels = function(meta, label_format) {
    matches = gregexpr("\\{[^{}]+\\}", label_format)
    fields = regmatches(label_format, matches)[[1]]
    parts = regmatches(label_format, matches, invert=TRUE)[[1]]
    columns = substring(fields, 2, nchar(fields)-1)
    missing = setdiff(columns, names(meta))
    if (length(missing) > 0)
        stop("Unknown sample label columns: ", paste(missing, collapse=", "))

    labels = rep(parts[1], nrow(meta))
    for (i in seq_along(columns))
        labels = paste0(labels, meta[[columns[i]]], parts[i+1], recycle0=TRUE)
    labels
}
