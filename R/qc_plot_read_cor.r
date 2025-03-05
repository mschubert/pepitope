calc_sample_cor = function(reps, keep=c("patient", "bc_type", "gene_name", "mut_id")) {
    merge_one = function(one_pat) {
        rdf1 = one_pat |> filter(rep == s_rep[1]) |>
            select(oligo_id, keep, smp1=smp, {{ r1 }}:=value)
        rdf2 = one_pat |> filter(rep == s_rep[2]) |>
            select(oligo_id, smp2=smp, {{ r2 }}:=value)
        re = inner_join(rdf1, rdf2, by="oligo_id", relationship="many-to-many") |>
            filter(! ({{ r1 }} == 0 & {{ r2 }} == 0))
    }
    s_rep = gtools::mixedsort(as.character(unique(reps$rep)))
    if (length(s_rep) == 1)
        return(list())
    r1 = paste("rep", s_rep[1])
    r2 = paste("rep", s_rep[2])
    re = lapply(split(reps, reps$patient), merge_one)
    attr(re, "reps") = c(r1, r2)
    re
}

plot_sample_cor = function(sample_cor_df, axes=c("rep 1", "rep 2"), add=list()) {
    ax1 = rlang::sym(axes[1])
    ax2 = rlang::sym(axes[2])
    if (is.data.frame(sample_cor_df)) {
        cors = sample_cor_df |> group_by(smp1, smp2) |>
            summarize(pearson = sprintf("R = %.3f", cor({{ ax1 }}, {{ ax2 }})))
        add = c(add, ggpp::geom_text_npc(data=cors, aes(label=pearson),
                                         npcx=0.1, npcy=0.9, size=2.5))
    }
    ggplot(sample_cor_df, aes(x={{ ax1 }}, y={{ ax2 }}, color=bc_type)) +
        geom_point(aes(text=gene_name), size=0.1) + add +
        scale_color_brewer(palette="Set1", guide=guide_legend(override.aes=list(size=2.5)), drop=FALSE) +
        scale_x_continuous(trans="log1p", breaks=c(0,20,500,10000,5e5), labels=c(0,20,500,"10k","500k")) +
        scale_y_continuous(trans="log1p", breaks=c(0,20,500,10000,5e5), labels=c(0,20,500,"10k","500k")) +
        facet_grid(rows=vars(smp2), cols=vars(smp1)) +
        ggtitle(paste("patient", unique(sample_cor_df$patient), "count correlations"))
}

plot_rep_cor = function(reps, pat) {
    if (length(unique(reps$rep)) != 3) stop("not implemented")
    cmp = reps |>
        filter(patient == pat) |>
        select(patient, origin, oligo_id, gene_name, bc_type, rep, value) |>
        tidyr::pivot_wider(names_from="rep", names_prefix="rep ", values_from="value", values_fill=0)
    df = list(cmp |> mutate(smp1 = "rep 1 vs 2", x = `rep 1`, y = `rep 2`),
              cmp |> mutate(smp1 = "rep 2 vs 3", x = `rep 2`, y = `rep 3`),
              cmp |> mutate(smp1 = "rep 1 vs 3", x = `rep 1`, y = `rep 3`)) |>
        bind_rows() |>
        select(smp1, smp2=origin, oligo_id, gene_name, bc_type, x, y)
    plot_sample_cor(df, axes=c("x", "y")) +
        coord_fixed() +
        labs(title = paste("patient", pat, "replicate correlations"),
             x = "replicates", y = "samples") +
        theme(strip.text.y = element_text(angle=0))
}
