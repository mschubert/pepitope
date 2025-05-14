#' Plot barcode overlap between different samples
#'
#' @param all_constructs  A named list of all constructs
#' @return  A ggplot2 object
#'
#' @export
plot_barcode_overlap = function(all_constructs, valid_barcodes) {
    dset = merge_constructs(all_constructs) |>
        select(bc_type, barcode) |>
        group_by(barcode) |>
            mutate(overlap = factor(n())) |>
        ungroup()

    if (!missing(valid_barcodes))
        dset = left_join(data.frame(barcode = valid_barcodes), dset)

    dset2 = dset |>
        mutate(barcode = as.integer(factor(barcode, levels=unique(barcode)))) |>
        filter(!is.na(bc_type))

    ggplot(dset2, aes(x=barcode, y=bc_type)) +
        geom_tile(aes(fill=overlap)) +
        scale_fill_brewer(palette="Set1") +
        ggtitle("Barcode overlap")
}
