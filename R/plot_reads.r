#' Plot the overall read counts
#'
#' @param dset  The `SummarizedExperiment` object from `count_bc`
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_layout
#' @export
plot_reads = function(dset) {
    plot_one = function(rsum, y, position) {
        ggplot(rsum, aes(x=label, y={{ y }}, fill=bc_type)) +
            geom_col(position=position) +
            scale_fill_brewer(palette="Set1", drop=FALSE) +
            coord_flip()
    }

    meta = as_tibble(colData(dset))
    read_summary = calc_representation(assay(dset), as_tibble(rowData(dset)), meta) |>
        group_by(sample_id, label, bc_type) |>
        summarize(reads=sum(value), nonzero_bcs=sum(value >= 10))
    unmapped = as_tibble(colData(dset)) |>
        mutate(bc_type=factor("unmapped"), reads=total_reads - mapped_reads) |>
        select(sample_id, label, bc_type, reads)
    both = bind_rows(read_summary, unmapped)

    noax = theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    plots = list(
        plot_one(both, reads, "stack") + ylab("Reads total") + ggtitle("Read counts"),
        plot_one(read_summary, reads, "fill") + noax + ylab("Mapped read fraction"),
        plot_one(read_summary, nonzero_bcs, "stack") +
            ylab("Barcodes total") + ggtitle("Barcodes >= 10 reads"),
        plot_one(read_summary, nonzero_bcs, "fill") + noax + ylab("Mapped barcode fraction")
    )
    wrap_plots(plots) + plot_layout(guides="collect") & theme(axis.title.y = element_blank())
}
