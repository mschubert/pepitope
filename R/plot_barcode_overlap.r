#' Plot barcode overlap between different samples
#'
#' @param all_constructs  A named list of all constructs
#' @return  A ggplot2 object
#' @export
plot_barcode_overlap = function(all_constructs, valid_barcodes) {
    dset = bind_rows(all_constructs, .id="bc_type") |>
        select(bc_type, barcode) |>
        group_by(barcode) |>
            mutate(overlap = factor(n())) |>
        ungroup()

    if (!missing(valid_barcodes))
        dset = left_join(data.frame(barcode = valid_barcodes), dset)

    dset$barcode = as.integer(factor(dset$barcode))
    ggplot(dset, aes(x=barcode, y=bc_type)) +
        geom_tile(aes(fill=overlap)) +
        scale_fill_brewer(palette="Set1") +
        ggtitle("Barcode overlap")
}
