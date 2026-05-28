#' Plot annotated read structure examples
#'
#' @param fq  Path to one FASTQ file
#' @param samples  A sample sheet as `data.frame` in tsv format. Requires the
#'      columns 'sample_id', 'patient', 'rep', 'origin', 'barcode'
#' @param all_constructs  A named list of all construct libraries
#' @param nrec  Number of FASTQ records to inspect
#' @return  A `ggplot2` object
#'
#' @import ggplot2
#' @importFrom BiocGenerics width
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

    structure = .rs_parse(ann$structure)
    features = rbind(
        data.frame(structure$sample, op=ifelse(structure$sample$revcomp, "B<", "B")),
        data.frame(structure$construct, op=ifelse(structure$construct$revcomp, "M<", "M"))
    )
    rects = merge(rows, features)

    ggplot() +
        geom_rect(data=rows, aes(xmin=0.5, xmax=read_width + 0.5, ymin=y - 0.2,
                                 ymax=y + 0.2),
                  fill="grey95", color="grey80") +
        geom_rect(data=rects, aes(xmin=start - 0.5, xmax=start + width - 0.5,
                                  ymin=y - 0.24, ymax=y + 0.24, fill=op),
                  alpha=0.45) +
        geom_text(data=bases, aes(x=position, y=y, label=base), size=3) +
        scale_y_reverse(breaks=rows$y, labels=paste0("read ", rows$read_idx),
                        expand=expansion(add=0.15)) +
        scale_x_continuous(breaks=seq.int(1L, read_width, by=max(1L, floor(read_width / 10L)))) +
        labs(x="Position", y=NULL, fill="Feature") +
        theme_minimal()
}
