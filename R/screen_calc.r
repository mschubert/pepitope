#' Calculate differential abundance of construct barcodes
#'
#' @param eset  A `SummarizedExperiment` object from `count_bc()`
#' @param comparisons  A character vector of sample and reference condition, or list thereof
#'
#' @importFrom stats setNames
#' @export
screen_calc = function(dset, comparisons) {
    if (length(unique(dset$patient)) > 1)
        stop("The 'patient' column can not span more than one value")

    eset = DESeq2::DESeqDataSet(dset, ~ rep + origin)

    eset$origin = factor(make.names(eset$origin))
    eset$rep = factor(eset$rep)
    DESeq2::sizeFactors(eset) = colSums(assay(dset)) / max(colSums(assay(dset)))
    mod = DESeq2::DESeq(eset)

    get_result = function(comp) {
        DESeq2::results(mod, contrast=c("origin", comp)) |>
            as.data.frame() |>
            tibble::rownames_to_column("barcode") |>
            as_tibble() |>
            filter(!is.na(log2FoldChange)) |>
            arrange(padj, pvalue) |>
            left_join(as.data.frame(rowData(eset)))
    }

    if (is.list(comparisons)) {
        lapply(comparisons, get_result) |>
            setNames(sapply(comparisons, paste, collapse=" vs "))
    } else {
        get_result(comparisons)
    }
}
