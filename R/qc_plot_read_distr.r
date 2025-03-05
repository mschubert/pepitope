calc_representation = function(lib_counts, bcs, meta) {
    reshape2::melt(lib_counts) |> as_tibble() |>
        dplyr::rename(oligo_id=Var1, sample_id=Var2) |>
        inner_join(bcs |> select(oligo_id, bc_type, gene_name, mut_id)) |>
        inner_join(meta) |>
        mutate(is_matched = bc_type == "TAA" |
               bc_type == as.character(sub("[^0-9]*([0-9]+).*", "\\1", patient))) |>
        group_by(sample_id) |>
            mutate(bcs_per_sample = n()) |>
        group_by(sample_id, bc_type, bcs_per_sample) |>
            mutate(frac = value / bcs_per_sample,
                   rnk = rank(value)/n()) |>
        ungroup()
}

plot_distr = function(reps) {
    ggplot(reps, aes(x=rnk, y=value, color=bc_type)) +
        geom_line(aes(linetype=is_matched)) +
        scale_linetype_manual(values=c("TRUE"="solid", "FALSE"="dashed")) +
        scale_color_brewer(palette="Set1", guide=guide_legend(override.aes=list(linewidth=0.8))) +
        scale_y_log10() +
        scale_x_continuous(limits=c(0,1), breaks=c(0,0.5,1)) +
        facet_wrap(~ short) +
        labs(x = "Ranked abundance",
             y = "Counts per barcode")
}
