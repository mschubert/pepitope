#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_layout
#' @export
plot_reads = function(reps, meta) {
    plot_one = function(rsum, y, position) {
        rsum$bc_type = factor(rsum$bc_type, levels=unique(both$bc_type))
        ggplot(rsum, aes(x=label, y={{ y }}, fill=bc_type)) +
            geom_col(position=position) +
            scale_fill_brewer(palette="Set1", drop=FALSE) +
            coord_flip()
    }

    reps_summary = reps |> group_by(sample_id, label, bc_type) |>
        summarize(reads=sum(value), nonzero_bcs=sum(value >= 10))
    unmapped = meta |> mutate(bc_type="unmapped", reads=total_reads - mapped_reads) |>
        select(sample_id, label, bc_type, reads)
    both = bind_rows(reps_summary, unmapped)

    noax = theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    plots = list(
        plot_one(both, reads, "stack") + ylab("Reads total") + ggtitle("Read counts") +
            geom_text(data=meta |> mutate(bc_type = NA), hjust=-0.1, size=2,
                      aes(y=total_reads, label=ifelse(mapped_reads < 1e5, "DISCARD", ""))),
        plot_one(reps_summary, reads, "fill") + noax + ylab("Mapped read fraction"),
        plot_one(reps_summary, nonzero_bcs, "stack") +
            ylab("Barcodes total") + ggtitle("Barcodes >= 10 reads"),
        plot_one(reps_summary, nonzero_bcs, "fill") + noax + ylab("Mapped barcode fraction")
    )
    wrap_plots(plots) + plot_layout(guides="collect") & theme(axis.title.y = element_blank())
}
