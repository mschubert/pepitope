#' Genomic track plot for variant + transcripts + RNA-seq
#'
#' @param rec A record
#' @param gene_id A gene id
#' @param res results object
#' @param gene gene name
#' @param txdb txdb object
#' @return A plot
#'
#' @importFrom GenomeInfoDb seqnames seqlevelsStyle
#' @importFrom Gviz GenomeAxisTrack AnnotationTrack GeneRegionTrack AlignmentsTrack plotTracks identifier<-
#' @importFrom IRanges start end
#' @keywords internal
plot_genomic_context = function(rec, gene_id, res, gene, txdb) {
#    gene_id = "ENSG00000162572"

    options(ucscChromosomeNames=FALSE) # needed by AlignmentsTrack
    seqlevelsStyle(txdb) = "ensembl"
    seqlevelsStyle(gene) = "ensembl"

    cur_var = res[res$GENEID == gene_id]
    cur_var$ref_nuc = NULL
    cur_var$alt_nuc = NULL
    cur_var$ref_prot = NULL
    cur_var$alt_prot = NULL
    cur_var = unique(cur_var)
    seqlevelsStyle(cur_var) = "ensembl"

    region = gene[gene$gene_id == gene_id]
    title = sprintf("%s (%s) @ chr%s:%i-%i", region$symbol, region$gene_id,
                    seqnames(region), start(region), end(region))
    tracks = list()
#    tracks$loc = IdeogramTrack(genome=genome(region)[1], chromosome=seqnames(region))
    tracks$axis = GenomeAxisTrack()
    tracks$var = AnnotationTrack(start=start(cur_var)-46, width=93, fill="red",
        chromosome=seqnames(cur_var), name="Variants",
        shape="ellipse", identifier=names(cur_var), showId=TRUE)
    identifier(tracks$var) = with(cur_var, sprintf("%s (p.%s%s%s) ",
        names(cur_var), REFAA, PROTEINLOC, VARAA))
    tracks$tx = GeneRegionTrack(txdb, name="Transcripts")
    tracks$rna = AlignmentsTrack(rec$til_rna$bam, type=c("coverage","sashimi"), name="RNA-seq")

    plotTracks(tracks, from=start(region), to=end(region), chromosome=seqnames(region),
               sizes=c(1,1,5,5), cex.main=1, main=title)
}
