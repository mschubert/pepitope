#' Plot barcode overlap between different samples
#'
#' @param all_constructs  A named list of all constructs
#' @param valid_barcodes  A character vector of possible barcodes (optional)
#' @return  A ggplot2 object
#'
#' @export
plot_barcode_overlap = function(all_constructs, valid_barcodes) {
    dset = merge_constructs(all_constructs, strict=FALSE) |>
        select(bc_type, barcode) |>
        group_by(barcode) |>
            mutate(overlap = factor(n())) |>
        ungroup()

    footer = list()
    if (!missing(valid_barcodes)) {
        n_invalid = length(setdiff(dset$barcode, valid_barcodes))
        dset = left_join(data.frame(barcode = valid_barcodes), dset)
        if (n_invalid == 0) {
            footer = list(labs(caption = "All barcodes valid"))
        } else {
            footer = list(labs(caption = paste(n_invalid, "invalid barcodes dropped")))
        }
    }

    dset2 = dset |>
        mutate(barcode = as.integer(factor(barcode, levels=unique(barcode)))) |>
        filter(!is.na(bc_type))

    ggplot(dset2, aes(x=barcode, y=bc_type)) +
        geom_tile(aes(fill=overlap)) +
        scale_fill_brewer(palette="Set1") +
        ggtitle("Barcode overlap") +
        footer
}
