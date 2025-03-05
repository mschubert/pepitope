library(dplyr)
library(ggplot2)
library(patchwork)

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

load_dset = function(args) {
    counts_stats = readr::read_tsv(args$stats_tsv) |> mutate(sample_id = sub("\\.R1$", "", label))
    counts_df = readr::read_tsv(args$counts_tsv)
    counts = data.matrix(counts_df[-(1:2)])
    rownames(counts) = counts_df$guide
    colnames(counts) = sub("\\.R1$", "", colnames(counts))
    meta = readr::read_tsv(args$meta_tsv) |>
        left_join(counts_stats |> select(sample_id, total_reads, mapped_reads)) |>
        mutate(patient = factor(patient),
               smp = paste(ifelse(nchar(origin)>8, stringr::word(origin, 1), origin), rep, sep="-"),
               short = paste(patient, smp),
               label = sprintf("%s (%s)", short, sample_id))
    list(meta=meta, counts=counts[,meta$sample_id, drop=FALSE])
}
