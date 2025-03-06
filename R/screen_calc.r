#' @importFrom stats setNames
calc_de = function(eset, cfg) {
    if (!is.null(cfg$patient))
        eset = eset[,eset$patient %in% cfg$patient]
    eset$origin = factor(make.names(sub("B cells", "Bcells", eset$origin)))
    eset$rep = factor(eset$rep)
    design(eset) = ~ rep + origin
    mod = DESeq2::estimateSizeFactors(eset) |> DESeq2::DESeq()

    get_result = function(comp) {
        DESeq2::results(mod, contrast=c("origin", comp)) |>
            as.data.frame() |>
            tibble::rownames_to_column("oligo_id") |>
            as_tibble() |>
            arrange(padj, pvalue) |>
            left_join(as.data.frame(SummarizedExperiment::rowData(eset))) |>
            mutate(bc_type = factor(bc_type))
    }
    lapply(cfg$comparisons, get_result) |>
        setNames(sapply(cfg$comparisons, paste, collapse=" vs "))
}

make_clean = function(res, sample, cap_fc=3) {
    res |>
        filter(bc_type %in% c(sub("+TAA", "", sample, fixed=TRUE), "shared", "TAA")) |>
        mutate(log2FoldChange = sign(log2FoldChange) * pmin(cap_fc, abs(log2FoldChange)))
}
