# make barcode lib
# https://github.com/hawkjo/freebarcodes
# https://www.pnas.org/doi/full/10.1073/pnas.1802640115
#' @importFrom tibble tibble
#' @importFrom utils read.table
make_lib = function() {
    url = "https://raw.githubusercontent.com/hawkjo/freebarcodes/master/barcodes/barcodes12-1.txt"
    bcs = read.table(url)$V1

    tsv = tibble::tibble(oligo_id = paste0("oligo_", seq_along(bcs)),
                 barcode = as.character(bcs)) |>
        dplyr::mutate(gene = oligo_id)
}
