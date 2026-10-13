#' Plot annotated read structure examples
#'
#' Forward and reverse barcode match counts are shown at each barcode start
#' position across all inspected reads, above up to ten read examples.
#' Only barcode types included in the inferred structure are annotated in the reads.
#' Read examples are shown without barcode annotations if the structure is unknown.
#'
#' @param fq  Path to one FASTQ file
#' @param samples  A sample sheet as `data.frame` in tsv format. Requires the
#'      columns 'sample_id', 'patient', 'rep', 'origin', 'barcode'
#' @param all_constructs  A named list of all construct libraries
#' @param nrec  Number of FASTQ records to inspect
#' @return  A patchwork object containing barcode counts and read examples
#'
#' @examples
#' if (interactive()) {
#'     samples = data.frame(sample_id="sample1", patient="pat1", rep="1",
#'         origin="library", barcode="GGG")
#'     constructs = list(pat1=data.frame(gene_name="GENE1", mut_id="GENE1_A1V",
#'         pep_id="GENE1_A1V", pep_type="alt", tiled="ATGGCCGCC", barcode_1="AAAA"))
#'     fq = tempfile(fileext=".fq")
#'     writeLines(c("@r1", "GGGAAAA", "+", "IIIIIII", "@r2", "GGGAAAA", "+", "IIIIIII"), fq)
#'     plot_read_structure(fq, samples, constructs, nrec=10)
#' }
#'
#' @import ggplot2
#' @importFrom BiocGenerics width
#' @importFrom patchwork wrap_plots
#' @importFrom tidyr pivot_longer
#' @export
plot_read_structure = function(fq, samples, all_constructs, nrec=100000L) {
    ann = .rs_annotate(fq, samples, all_constructs, nrec)
    n_reads = min(10L, length(ann$reads))
    read_idx = sample(seq_along(ann$reads), n_reads)
    reads = ann$reads[read_idx]
    read_width = width(reads)[1]
    rows = data.frame(read=seq_len(n_reads), y=seq_len(n_reads) * 0.55, read_idx=read_idx,
                      read_width=read_width)
    bases = do.call(rbind, Map(function(read, seq) {
        data.frame(read=read, position=seq_len(read_width), base=strsplit(seq, "")[[1]])
    }, rows$read, as.character(reads)))
    bases = merge(bases, rows[c("read", "y")])

    known_structure = ann$structure != "unknown"
    feature_labels = c(B="Sample (B)", M="Construct (M)")
    position_scale = scale_x_continuous(limits=c(0.5, read_width + 0.5),
        breaks=seq.int(1L, read_width, by=max(1L, floor(read_width / 10L))),
        expand=expansion(mult=0))
    counts = ann$counts
    counts$position = seq_len(nrow(counts))
    counts = pivot_longer(counts, cols=c("B", "^B", "M", "^M"), names_to="op", values_to="count")
    counts$type = sub("^\\^", "", counts$op)
    counts$orientation = factor(ifelse(startsWith(counts$op, "^"), "Reverse", "Forward"),
                                levels=c("Forward", "Reverse"))
    count_plot = ggplot(counts, aes(x=position, y=count, colour=type)) +
        geom_line(linewidth=0.4) +
        geom_point(size=0.9) +
        facet_grid(orientation ~ .) +
        scale_colour_discrete(name="Feature", limits=names(feature_labels), labels=feature_labels,
                              guide=if (known_structure) "none" else "legend") +
        scale_y_continuous(limits=c(0, NA)) +
        position_scale +
        labs(x=NULL, y="Barcode matches",
             title=paste("Read structure identified as:", ann$structure),
             subtitle=paste("First", format(nrec, big.mark=","), "reads inspected")) +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

    matches = .rs_match(reads, ann$structure, .read_samples(samples)$barcode,
                        merge_constructs(all_constructs)$barcode)
    barcode_layers = NULL
    if (known_structure) {
        rects = merge(rows, matches, by="read")
        rects$barcode = ifelse(rects$present, "Present", "Absent")
        rects$fill = ifelse(rects$present, rects$type, NA_character_)
        barcode_layers = list(
            geom_rect(data=rects, aes(xmin=start - 0.5, xmax=start + width - 0.5,
                                     ymin=y - 0.24, ymax=y + 0.24, fill=fill,
                                     colour=type, linetype=barcode), alpha=0.45, show.legend=TRUE),
            scale_colour_discrete(aesthetics=c("colour", "fill"), name="Feature",
                                  limits=names(feature_labels), labels=feature_labels,
                                  na.value=NA, na.translate=FALSE),
            scale_linetype_manual(name="Barcode", values=c(Present="solid", Absent="solid"),
                                  limits=c("Present", "Absent"),
                                  guide=guide_legend(override.aes=list(fill=c("grey40", NA),
                                                                      colour="grey40")))
        )
    }

    read_plot = ggplot() +
        geom_rect(data=rows, aes(xmin=0.5, xmax=read_width + 0.5, ymin=y - 0.2, ymax=y + 0.2),
                  fill="grey95", color="grey80") +
        barcode_layers +
        geom_text(data=bases, aes(x=position, y=y, label=base), size=3) +
        scale_y_reverse(breaks=rows$y, labels=paste0("read ", rows$read_idx),
                        expand=expansion(add=0.15)) +
        position_scale +
        labs(x="Position", y=NULL)

    wrap_plots(count_plot, read_plot, ncol=1, heights=c(2, 4), guides="collect") & theme_minimal()
}
