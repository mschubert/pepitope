context("annotation")

fname = system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vr = VariantAnnotation::readVcfAsVRanges(fname)
ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38

is_prot = GeneBiotypeFilter("protein_coding")
is_chr = SeqNameFilter(c(1:22,'X','Y'))
tx = transcripts(ens106, filter=c(is_prot, is_chr))
tx_coding = names(tx)[tx$tx_biotype=="protein_coding"]

#res = predictCoding(vr, ens106, asm)
res = suppressWarnings(annotate_coding(vr, ens106, asm, tx_coding))
