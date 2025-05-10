#' Plot screen results
#'
#' @param res     A results `data.frame` from `screen_calc()`
#' @param sample  Which library to plot (default: all)
#' @param links   Whether to draw arrows between ref and significant alt peptides
#' @param labs    Whether to label genes in less dense areas
#' @param cap_fc  Maximum amount of fold-change to limit values to
#' @return  A `ggplot2` object of the differential expression results
#'
#' @export
plot_screen = function(res, sample, links=TRUE, labs=TRUE, cap_fc=8) {
    res$log2FoldChange = sign(res$log2FoldChange) * pmin(abs(res$log2FoldChange), cap_fc)
    lab = res |> filter(padj<0.1) |> group_by(gene_name) |>
        filter(n_distinct(pep_type)==1 | is.na(pep_type) | pep_type=="alt")

    p = ggplot(res, aes(x=baseMean, y=log2FoldChange)) +
        geom_hline(yintercept=0, color="grey") +
        geom_point(aes(color=bc_type, shape=pep_type, size=padj<0.1, alpha=padj<0.1)) +
        scale_shape_manual(values=c(ref=1, alt=19), na.value=19) +
        scale_size_manual(values=c("TRUE"=1, "FALSE"=0.2), na.value=0.2) +
        scale_alpha_manual(values=c("TRUE"=0.5, "FALSE"=0.2), na.value=0.2) +
        scale_color_brewer(palette="Set1", drop=FALSE) +
        scale_x_log10(limits=c(10,NA))
    if (links) {
        link = make_links(res, sample)
        if (nrow(link) > 0)
            p = p + ggpp::stat_dens2d_filter(data=link, geom="segment", color="black", alpha=0.2, linewidth=0.3,
                aes(x=baseMean_alt, xend=baseMean_ref, y=log2FoldChange_alt, yend=log2FoldChange_ref),
                arrow=arrow(length=unit(0.01, "npc"), type="closed", ends="first"),
                keep.fraction=1, keep.number=30)
    }
    if (labs && nrow(lab) > 0) {
        p = p + ggpp::stat_dens2d_filter_g(data=lab, aes(label=gene_name, color=bc_type),
            geom="text_repel", keep.fraction=1, keep.number=45,
            size=1.5, min.segment.length=0, segment.alpha=0.3, segment.size=0.3)
    }
    p
}

#' @keywords internal
make_links = function(res, sample) {
    if (missing(sample))
        sample = unique(res$bc_type)
    arrs = res |> filter(bc_type %in% sample) |>
        select(mut_id, pep_type, baseMean, log2FoldChange, padj)
    ar1 = arrs |> filter(pep_type == "ref") |> select(-pep_type) |>
        dplyr::rename(baseMean_ref=baseMean, log2FoldChange_ref=log2FoldChange)
    ar2 = arrs |> filter(pep_type == "alt", padj<0.1) |> select(-pep_type) |>
        dplyr::rename(baseMean_alt=baseMean, log2FoldChange_alt=log2FoldChange, padj_alt=padj)
    inner_join(ar1, ar2, relationship="many-to-many")
}
