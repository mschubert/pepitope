#' Plot the overall read counts
#'
#' Reads with a recognized construct barcode but an unmatched sample barcode are
#' included as an "Unknown sample" group using `rowData(dset)$unassigned_reads`.
#' Reads matching neither barcode are included as unmapped reads in this group using
#' `metadata(dset)$unmapped_reads`.
#'
#' @param dset  The `SummarizedExperiment` object from `count_fastq`
#' @param n_cutoff  Minimum reads for counting a barcode as represented
#' @param label_format  Sample label template with column placeholders
#' @return  A patchwork object containing the read count summary plots
#'
#' @examples
#' if (interactive()) {
#'     samples = data.frame(sample_id="sample1", minigenes="pat1", rep="1",
#'         origin="library", barcode="GGG")
#'     constructs = list(pat1=data.frame(gene_name="GENE1", mut_id="GENE1_A1V",
#'         pep_id="GENE1_A1V", pep_type="alt", tiled="ATGGCCGCC", barcode_1="AAAA"))
#'     fq = tempfile(fileext=".fq")
#'     writeLines(c("@r1", "GGGAAAA", "+", "IIIIIII", "@r2", "GGGAAAA", "+", "IIIIIII"), fq)
#'     dset = count_fastq(fq, samples, constructs, read_structure="3B4M", verbose=FALSE)
#'     plot_read_count(dset)
#' }
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom S4Vectors metadata
#' @export
plot_read_count = function(dset, n_cutoff=10, label_format="{sample_id} ({minigenes})") {
    plot_one = function(rsum, y, position) {
        ggplot(rsum, aes(x=sample_id, y={{ y }}, fill=minigenes)) +
            geom_col(position=position, show.legend=c(unused=TRUE, unmapped=TRUE)) +
            scale_fill_manual(values=make_pal(rsum$minigenes), drop=FALSE) +
            scale_x_discrete(limits=rev(meta$sample_id),
                              labels=setNames(meta$label, meta$sample_id)) +
            coord_flip()
    }

    bcs = as_tibble(rowData(dset))
    counts = cbind(assay(dset), unassigned=bcs$unassigned_reads)
    meta = as_tibble(colData(dset))
    meta$label = .format_sample_labels(meta, label_format)
    add_meta = tibble(sample_id="unassigned", label="Unknown sample(s)",
                      total_reads=sum(bcs$unassigned_reads) + metadata(dset)$unmapped_reads,
                      mapped_reads=sum(bcs$unassigned_reads))
    meta = bind_rows(meta, add_meta)
    read_summary = calc_representation(counts, bcs, meta) |>
        group_by(sample_id, label, minigenes) |>
        summarize(reads=sum(value), nonzero_bcs=sum(value >= n_cutoff))
    unmapped = meta |>
        mutate(minigenes=factor("unmapped"), reads=total_reads - mapped_reads) |>
        select(sample_id, label, minigenes, reads)
    both = bind_rows(read_summary, unmapped)
    levels(read_summary$minigenes) = levels(both$minigenes)

    noax = theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    plots = list(
        plot_one(both, reads, "stack") + ylab("Reads total") + ggtitle("Read counts"),
        plot_one(read_summary, reads, "fill") + noax + ylab("Mapped read fraction"),
        plot_one(read_summary, nonzero_bcs, "stack") +
            ylab("Barcodes total") + ggtitle(sprintf("Barcodes >= %s reads", n_cutoff)),
        plot_one(read_summary, nonzero_bcs, "fill") + noax + ylab("Mapped barcode fraction")
    )
    wrap_plots(plots) + plot_layout(guides="collect") & theme(axis.title.y = element_blank())
}

#' @rdname plot_read_count
#' @export
plot_reads = function(dset) {
    warning("This function will be removed. Use plot_read_count() instead")
    plot_read_count(dset)
}
