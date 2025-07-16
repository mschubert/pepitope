#' Bind the tables of all library annotations together
#'
#' @param all_constructs  A named list of all construct libraries
#' @param strict  Whether to make sure the merged table is valid
#' @return  A `data.frame` with all constructs
#'
#' @keywords internal
merge_constructs = function(all_constructs, strict=TRUE) {
    if (!is.list(all_constructs) || is.data.frame(all_constructs))
        stop("'all_constructs' needs to be a named list but is not a list")
    if (is.null(names(all_constructs)))
        stop("'all_constructs' needs to be a named list but has no names")

    expand_barcode = function(name, df) {
        if ("mut_id" %in% colnames(df) && ! "gene_name" %in% colnames(df))
            df$gene_name = sub("_[^_]+$", "", df$mut_id)
        if ("gene_name" %in% colnames(df) && ! "mut_id" %in% colnames(df))
            df$mut_id = df$gene_name
        if ("type" %in% colnames(df) && ! "pep_type" %in% colnames(df))
            df = dplyr::rename(df, pep_type=type)
        req = c("gene_name", "mut_id", "pep_id", "tiled") # pep_type can be NA
        bc_fields = grep("^barcode", colnames(df), value=TRUE)
        if (all(req %in% colnames(df)) && length(bc_fields) > 0) {
            tidyr::pivot_longer(df, bc_fields, values_to="barcode") |> select(-name)
        } else {
            missing = paste(setdiff(req, colnames(df)), collapse=", ")
            if (length(bc_fields) == 0)
                missing = c(missing, "barcode")
            warning("Skipping ", sQuote(name), " because of missing fields: ", missing)
        }
    }
    construct_df = mapply(expand_barcode, name=names(all_constructs),
                          df=all_constructs, SIMPLIFY=FALSE) |>
        bind_rows(.id="bc_type")

    if (strict) {
        if (nrow(construct_df) == 0)
            stop("No data left after merging 'all_constructs'")
        if (any(duplicated(construct_df$barcode)))
            stop("Duplicate barcodes are not allowed but found in 'all_constructs'")
    }

    construct_df
}
