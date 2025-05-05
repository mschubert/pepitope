#' Calculate ranked barcode representation across samples
#'
#' @param lib_counts  A count matrix
#' @param bcs  Barcodes
#' @param meta  Sample sheet and read count summaries
#' @keywords internal
calc_representation = function(lib_counts, bcs, meta) {
    reshape2::melt(lib_counts) |> as_tibble() |>
        dplyr::rename(barcode=Var1, sample_id=Var2) |>
        inner_join(bcs |> select(barcode, bc_type, gene_name, mut_id, pep_id, pep_type)) |>
        inner_join(meta |> select(-barcode)) |>
        rowwise() |>
            mutate(is_matched = bc_type %in% strsplit(as.character(patient), "+", fixed=TRUE)[[1]]) |>
        group_by(sample_id) |>
            mutate(bcs_per_sample = n()) |>
        group_by(sample_id, bc_type, bcs_per_sample) |>
            mutate(frac = value / bcs_per_sample,
                   rnk = (rank(value, ties.method="first")-1) / (n()-1),
                   rlab = paste0(rank(value, ties.method="first"), "/", n())) |>
        ungroup()
}

#' Plot the read distribution across barcodes
#'
#' @param dset  The `SummarizedExperiment` object from `count_bc`
#' @return  A `ggplot2` object
#' @export
plot_distr = function(dset) {
    reps = calc_representation(assay(dset), as_tibble(rowData(dset)), as_tibble(colData(dset))) |>
        filter(is_matched | value > 0) |>
        mutate(text = sprintf("%s (%s)\n%i reads (%.1f%%)\n%s (%s)",
            ifelse(is.na(pep_id), oligo_id, pep_id),
            ifelse(is.na(pep_type), "unused", pep_type),
            as.integer(value), frac,
            ifelse(is_matched, "matched", "contamination"), rlab))

    ggplot(reps, aes(x=rnk, y=value, group=bc_type, color=bc_type, text=text)) +
        geom_line(aes(linetype=is_matched)) +
        scale_linetype_manual(values=c("TRUE"="solid", "FALSE"="dashed")) +
        scale_color_manual(values=make_pal(reps$bc_type),
                           guide=guide_legend(override.aes=list(linewidth=0.8))) +
        scale_y_continuous(trans="log1p", breaks=c("0"=0,"100"=100,"10k"=1e4,"1M"=1e6)) +
        scale_x_continuous(limits=c(0,1), breaks=c(0,0.5,1)) +
        facet_wrap(~ short) +
        labs(x = "Ranked abundance",
             y = "Counts per barcode")
}

#' @keywords internal
make_pal = function(values) {
    pats = setdiff(na.omit(values), c("unused", "unmapped"))
    if (length(pats) <= 9) {
        cols = RColorBrewer::brewer.pal(length(pats), "Set1")
        col_pats = setNames(cols, pats)
    } else {
        col_fun = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
        cols = col_fun(length(pats)+1)
        col_pats = setNames(cols[-length(cols)], pats)
    }
    c(col_pats, unused="grey80", unmapped="grey30")
}
