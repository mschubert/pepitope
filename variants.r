import_package("GenomicFeatures", attach=TRUE)
import_package("VariantAnnotation", attach=TRUE)
import_package("dplyr", attach=TRUE)
sys = modules::import('sys')
subseq = Biostrings::subseq
nchar = Biostrings::nchar
vcountPattern = Biostrings::vcountPattern
vmatchPattern = Biostrings::vmatchPattern

#' Annotate VCF variants with coding changes
#'
#' @param rec   A sample record object
#' @param txdb  Txdb object
#' @param asm   Genomic sequence BSGenome object
#' @param tx_coding  Character vector of ENST0000 IDs that are protein coding
#' @param tumor_cov  Part of column name to get ref/alt coverage (regex)
#' @return      A GRanges object with annotated variants
annotate_coding = function(rec, txdb, asm, tx_coding, tumor_cov="tumor_DNA") {
    vr = readVcfAsVRanges(rec$dna$vcf_diff, "GRCh38")
    vr = vr[grepl(tumor_cov, sampleNames(vr))]
    vr$sampleNames = sampleNames(vr)
    vr$ref = ref(vr)
    vr$alt = alt(vr)
    vr$cov_ref = refDepth(vr)
    vr$cov_alt = altDepth(vr)
    codv = predictCoding(vr, txdb, asm)

#    codv2 = predictCoding(vcf_tumor, txdb, asm) # 1637 var names @codv, 26k @codv2, 223 common?!
#    splice = locateVariants(vcf_diff, txdb, SpliceSiteVariants())

    tx = transcripts(txdb)
    utr3 = threeUTRsByTranscript(txdb)
    codv$tx_name = tx$tx_name[match(codv$TXID, tx$tx_id)]
    codv = codv[codv$tx_name %in% tx_coding]
    coding_ranges = cdsBy(txdb)[codv$TXID]
    codv$ref_nuc = extractTranscriptSeqs(asm, coding_ranges)
    codv$ref_prot = Biostrings::translate(codv$ref_nuc)

    # filter for proper ORFs
    has_start = subseq(codv$ref_prot,1,1) == "M"
    has_stop = subseq(IRanges::reverse(codv$ref_prot),1,1) == "*"
    n_stop = vcountPattern("*", codv$ref_prot)
    codv = codv[has_start & has_stop & n_stop==1] #TODO: replace alt init codons by M

    # get coding sequences with updated variants
    upd_seqs = function(i) {
        new = codv$varAllele[[i]]
        Biostrings::replaceAt(codv$ref_nuc[[i]], codv$CDSLOC[i], new)
    }
#    atn = Biostrings::replaceAt(codv$ref_nuc, as(codv$CDSLOC, "IRangesList"), codv$varAllele)
    codv$alt_nuc = Biostrings::DNAStringSet(lapply(seq_along(codv$ref_nuc), upd_seqs))

    # get protein sequences and adjust nuc for premature stop
    codv$alt_prot = Biostrings::translate(codv$alt_nuc)
    stops = vmatchPattern("*", codv$alt_prot)
    first = sapply(stops, function(s) IRanges::start(s)[1])
    changed = which(is.na(first) | first != nchar(codv$alt_prot))
    for (i in changed) {
        if (is.na(first[i])) {
            if (! codv$TXID[[i]] %in% names(utr3))
                next
            nuc_utr3 = getSeq(asm, utr3[[codv$TXID[i]]])
            codv$alt_nuc[i] = Biostrings::xscat(codv$alt_nuc[i], nuc_utr3)
            codv$alt_prot[i] = Biostrings::translate(codv$alt_nuc[i])
            first[i] = IRanges::start(vmatchPattern("*", codv$alt_prot[i])[[1]])[1]
        }
        codv$alt_prot[i] = subseq(codv$alt_prot[i], 1, first[i])
        codv$alt_nuc[i] = subseq(codv$alt_nuc[i], 1, first[i]*3)
    }

    check_silent = function(ref, alt) {
        offsets = rep(0, length(ref))
        chk = seq_along(ref)
        for (i in seq_len(max(nchar(ref)))) {
            chk = chk[nchar(ref[chk]) >= i & nchar(alt[chk]) >= i]
            s_mtch = substr(ref[chk],i,i) == substr(alt[chk],i,i)
            if (! any(s_mtch))
                break
            chk = chk[s_mtch]
            offsets[chk] = offsets[chk] + 1
        }
        offsets
    }
    codv$silent_start = check_silent(codv$REFAA, codv$VARAA)
    codv$silent_end = check_silent(Biostrings::reverse(codv$REFAA), Biostrings::reverse(codv$VARAA))
    codv$silent_end[codv$silent_start == nchar(codv$REFAA) | # eg. A>ASA not double
                    codv$silent_start == nchar(codv$VARAA)] = 0 # eg. AA>A not double
    silent = codv$silent_start + codv$silent_end

    codv$CONSEQUENCE = as.character(codv$CONSEQUENCE)
    codv$CONSEQUENCE[IRanges::start(codv$CDSLOC) == 1 & codv$VARCODON != "ATG"] = "nostart"
    codv$CONSEQUENCE[silent == nchar(codv$REFAA) & nchar(codv$VARAA) > nchar(codv$REFAA)] = "insertion"
    codv$CONSEQUENCE[silent == nchar(codv$VARAA) & nchar(codv$REFAA) > nchar(codv$VARAA)] = "deletion"

    # check if we didn't change the length of any nuc
    stopifnot(with(codv,
        nchar(codv$ref_nuc) - nchar(codv$REFCODON) == nchar(codv$alt_nuc) - nchar(codv$VARCODON) |
        codv$CONSEQUENCE %in% c("frameshift", "nonsense", "nostart") |
        vcountPattern("*", codv$REFAA) > 0 | # transcript extension
        vcountPattern("*", codv$VARAA) > 0 # "missense"+nonsense
    ))
    codv
}

#' Subset nucleotide/protein sequences to codon +/- 45 bp context
#'
#' @param codv  Annotated variants from `annotate_coding()`
#' @param ctx_codons  How many flanking codons each to include in the context
#' @return  GRanges object with sequence information of only context
subset_context = function(codv, ctx_codons=15) {
    ctx = ctx_codons * 3
    stopAA = sapply(vmatchPattern("*", codv$VARAA), function(x) IRanges::start(x)[1]-1) * 3
    len_delta = pmin(nchar(codv$VARCODON), stopAA, na.rm=TRUE) - nchar(codv$REFCODON)

    roi_codon_start = floor((IRanges::start(codv$CDSLOC)-1)/3 + codv$silent_start) * 3 + 1
    roi_codon_end_ref = ceiling(IRanges::end(codv$CDSLOC)/3 - codv$silent_end) * 3
    len_ref = nchar(codv$ref_nuc) - 3
    len_alt = nchar(codv$alt_nuc) - 3

    is_frameshift = abs(len_delta) %% 3 != 0
    roi_codon_end_alt = pmin(ceiling((IRanges::end(codv$CDSLOC)+len_delta)/3 - codv$silent_end) * 3, len_alt)
    roi_codon_end_alt[is_frameshift] = nchar(codv$alt_nuc)[is_frameshift] - 3

    ctx_start = pmax(1, roi_codon_start - ctx)
    ctx_end_ref = pmin(len_ref, roi_codon_end_ref + ctx)
    ctx_end_alt = pmin(floor(len_alt/3)*3, roi_codon_end_alt + ctx) #TODO: better to force %%3 for nuc length?
    ctx_len_alt = ctx_end_alt - ctx_start + 1

    ctx_start_over = pmax(0, ctx + 1 - roi_codon_start)
    ctx_end_over = pmax(0, roi_codon_end_alt + ctx - len_alt)

    ext_by = pmin(len_alt, 2*ctx+3) - ctx_len_alt
    add_to_end = pmin(ext_by, ctx_start_over) %>% pmin(len_alt - ctx_start) %>% pmax(0)
    add_to_start = pmin(ext_by, ctx_end_over) %>% pmin(ctx_start-1) %>% pmax(0)
    add_to_start[codv$VARAA == "*"] = 0
    stopifnot(add_to_start == 0 | add_to_end == 0,
              c(add_to_start, add_to_end) %% 3 == 0,
              c(add_to_start, add_to_end) <= ctx+3)

    # get nuc and protein sequences incl context
    codv$ref_nuc = subseq(codv$ref_nuc, ctx_start-add_to_start, ctx_end_ref+add_to_end) #TODO: merge if context overlaps? [only correct if same allele/read support]
    codv$alt_nuc = subseq(codv$alt_nuc, ctx_start-add_to_start, ctx_end_alt+add_to_end)
    codv$ref_prot = Biostrings::translate(codv$ref_nuc) #TODO: no alt init codons (atn M replaced by ATG above)
    codv$alt_prot = Biostrings::translate(codv$alt_nuc)
    codv$alt_nnuc = nchar(codv$alt_nuc)
    codv$alt_shift = add_to_end - add_to_start

    codv
}

#' Add gene expression values to variant result
#'
#' @param res  Annotated variants from `annotate_coding()`
#' @param rec  Sample record from config file
#' @return     Annotated variants including gene counts and TPM
add_gex = function(res, rec) {
    counts = readr::read_tsv(rec$til_rna$count, col_names=c("gene_id", "gene_name", "count"))
    gex = readr::read_tsv(rec$til_rna$tpm) %>%
        inner_join(counts) %>%
        dplyr::select(gene_id, gene_name, locus, count, TPM)

    res$gene_count = gex$count[match(res$GENEID, gex$gene_id)]
    res$gene_tpm = gex$TPM[match(res$GENEID, gex$gene_id)]
    res
}

#' Genomic track plot for variant + transcripts + RNA-seq
plot_genomic_context = function(rec, gene_id, res, gene, txdb) {
#    gene_id = "ENSG00000162572"

    library(Gviz)
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

#' Save results as xlsx sheets (full, filtered, peptides)
#'
#' @param res  A full results GRanges object from `annotate_coding()`
#' @param fname  File name to save results to
save_xlsx = function(res, fname, min_cov=2, min_af=0.1) {
    gr2df = function(gr) as_tibble(as.data.frame(unname(gr))) %>%
        mutate(var_id = names(gr)) %>%
        dplyr::select(var_id, everything()) %>%
        arrange(order(gtools::mixedorder(var_id)))

    # changes peptide, is unique and is expressed
    names(res) = sprintf("%s:%i", GenomicRanges::seqnames(res), IRanges::start(res))
    subs = subset_context(res[! res$CONSEQUENCE %in% c("synonymous", "nonsense", "nostart")])
    alt_in_ref = function(a,r) grepl(as.character(a), as.character(r), fixed=TRUE)
    subs = subs[!mapply(alt_in_ref, a=subs$alt_prot, r=subs$ref_prot)]
    if ("gene_count" %in% colnames(S4Vectors::mcols(subs)))
        subs = subs[!is.na(subs$gene_count) & subs$gene_count > 0] # & subs$gene_tpm > 0]
    subs = subs[subs$cov_alt/subs$cov_ref >= min_af & subs$cov_alt >= min_cov]
    subs = subs[!duplicated(subs$alt_prot)]

    # peptide is not contained within another peptide
    contained_in = function(i) any(grepl(ac2[i], ac2[-i], fixed=TRUE))
    ac = as.character(subs$alt_prot)
    ac2 = substr(ac, pmax(1, 1-subs$alt_shift/3), nchar(ac)-pmax(0, subs$alt_shift/3))
    any_con = sapply(seq_along(ac2), contained_in)
    subs = subs[!any_con]

    table(nchar(subs$alt_nuc))
    stopifnot(nchar(subs$alt_nuc) %% 3 == 0)

    # tile peptides to have max 93 nt length
    tile_cDNA = function(p) {
        ntile = ceiling(nchar(p)/93)
        starts = round(seq(0, nchar(p)-93, length.out=ntile)/3) * 3 + 1
        lapply(starts, function(s) substr(p, s, s+92))
    }
    pep = with(subs, tibble(var_id=names(subs), gene_id=subs$GENEID,
                            tx_id=subs$tx_name, ref=as.character(subs$ref_nuc),
                            alt=as.character(subs$alt_nuc))) %>%
        tidyr::pivot_longer(c(ref, alt), names_to="type", values_to="cDNA") %>%
        rowwise() %>% mutate(tiled = list(tile_cDNA(cDNA))) %>% ungroup() %>%
        mutate(n_tiles = sapply(tiled, length)) %>%
        tidyr::unnest(tiled) %>%
        mutate(tiled = unlist(tiled, use.names=FALSE),
               nt = nchar(tiled),
               peptide = as.character(Biostrings::translate(Biostrings::DNAStringSet(tiled), no.init.codon=TRUE)))

#    stopifnot(pep$type[duplicated(pep$tiled)] == "ref")
    pep = pep[!duplicated(pep$tiled),]

    pep$BbsI_replaced = vcountPattern("GAAGAC", pep$tiled) + vcountPattern("GTCTTC", pep$tiled)
    pep$tiled = sapply(pep$tiled, remove_cutsite, site="GAAGAC", seed=178529, USE.NAMES=FALSE)
#    stopifnot(pep$peptide == as.character(Biostrings::translate(Biostrings::DNAStringSet(pep$tiled), no.init.codon=TRUE)))

    sv = list(
        `All Variants` = res %>% select(-CDSLOC) %>% gr2df() %>%
            dplyr::select(-tx_name, -(ref_nuc:alt_prot)) %>% distinct(),
        `Unique Protein-Coding` = gr2df(subs),
        `93 nt Peptides` = pep %>% dplyr::select(var_id:cDNA, n_tiles, BbsI_replaced, tiled, nt, peptide)
    )
    writexl::write_xlsx(sv, path=fname)
}

#' Remove a Restriction Enzyme cut site but keep AA
#'
#' @param nuc   cDNA nucleotide string
#' @param site  Recognition site to be replaced (fwd+rev comp)
#' @param seed  Set random seed to select same changes on multiple runs
#' @return      cDNA with minimal changes to no longer contain the cut site
remove_cutsite = function(nuc, site, seed=NULL) {
    set.seed(seed)
    revtrans = function(aa) {
        aa_split = strsplit(as.character(aa), "+")[[1]]
        tr = Biostrings::getGeneticCode()
        lapply(aa_split, function(x) names(tr)[tr == x]) %>%
            do.call(tidyr::crossing, .) %>%
            rowwise() %>%
            purrr::pmap_chr(paste0)
    }
    alt_nuc = function(nuc, match) {
        subs = subseq(Biostrings::DNAString(nuc), IRanges::start(match), IRanges::end(match))
        possib = revtrans(Biostrings::translate(subs, no.init.codon=TRUE))
        lapply(possib, Biostrings::replaceAt, x=nuc, at=match)
    }

    rc = function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
    m = c(vmatchPattern(site, nuc), vmatchPattern(rc(site), nuc)) %>% unlist()
    if (length(m) == 0)
        return(nuc)

    IRanges::start(m) = floor((IRanges::start(m)-1)/3) * 3 + 1
    IRanges::end(m) = ceiling((IRanges::end(m))/3) * 3

    nucs = nuc = Biostrings::DNAStringSet(nuc)
    for (i in seq_along(m))
        nucs = lapply(nucs, alt_nuc, match=m[i]) %>% do.call(c, .)
    nucs = Biostrings::DNAStringSet(nucs)

    valid = vcountPattern(site, nucs) + vcountPattern(rc(site), nucs) == 0
    nchange = mapply(function(i) stringdist::stringdist(subseq(nuc, m[i]), subseq(nucs, m[i])),
                     i=seq_along(m)) %>% as.matrix() %>% rowSums()
    min_valid = which(valid & nchange == min(nchange[valid]))
    as.character(nucs[sample(min_valid, 1)])
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', 'samples.yaml'),
        opt('a', 'patient', 'identifier', 'M14TIL102'),
        opt('o', 'outfile', 'rds', 'M14TIL102.rds'),
        opt('x', 'export', 'xlsx', 'M14TIL102.xlsx')
#        opt('p', 'plotfile', 'pdf', 'var_M14TIL102.pdf')
    )

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
})
