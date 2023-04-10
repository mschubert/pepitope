#' Main function, remove this later
main = function(args) {
#    sys = modules::import('sys')
#    args = sys$cmd$parse(
#        opt('c', 'config', 'yaml', 'samples.yaml'),
#        opt('a', 'patient', 'identifier', 'M14TIL102'),
#        opt('o', 'outfile', 'rds', 'M14TIL102.rds'),
#        opt('x', 'export', 'xlsx', 'M14TIL102.xlsx')
##        opt('p', 'plotfile', 'pdf', 'var_M14TIL102.pdf')
#    )

    cfg = yaml::read_yaml(args$config)
    rec = cfg$samples[[args$patient]]
    rec$dna = lapply(rec$dna, function(p) file.path(cfg$path$dna, p))
    rec$til_rna = lapply(rec$til_rna, function(p) file.path(cfg$path$til_rna, p))

    # AnnotationHub::query(AnnotationHub::AnnotationHub(), c("EnsDb", "sapiens"))
    ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
    seqlevelsStyle(ens106) = "UCSC"
    is_prot = GeneBiotypeFilter("protein_coding")
    is_chr = SeqNameFilter(c(1:22,'X','Y'))
    gene = genes(ens106, filter=c(is_prot, is_chr))
    tx = transcripts(ens106, filter=c(is_prot, is_chr))
    tx_coding = names(tx)[tx$tx_biotype=="protein_coding"]

    asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    seqlevelsStyle(asm) = "UCSC"
    txdb = makeTxDbFromGFF("data/Homo_sapiens.GRCh38.106.chr.gtf.gz")
    seqlevelsStyle(txdb) = "UCSC"

    res = annotate_coding(rec, txdb, asm, tx_coding) %>%
        add_gex(rec)

    saveRDS(res, file=args$outfile)
    save_xlsx(res, args$export)
#todo: use canonical transcript if all have same ctx

#    cur_var2 = res3[res3$GENEID =="ENSG00000162572"]
#    algn = pairwiseAlignment(cur_var2$ref_nuc, cur_var2$alt_nuc)
}
