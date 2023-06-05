#' Save results as xlsx sheets (full, filtered, peptides)
#'
#' @param res      A full results GRanges object from `annotate_coding()`
#' @param fname    File name to save results to
#' @param min_cov  Minimum number of reads to span the ALT allele
#' @param min_af   Minimum allele frequency of the ALT allele
#' @param tile_size  Oligo tiling size
#' @param tile_ov    Oligo tiling overlap
#'
#' @importFrom dplyr `%>%` rowwise mutate select arrange group_by ungroup as_tibble
#' @importFrom Biostrings nchar translate DNAStringSet vcountPattern
save_xlsx = function(res, fname, min_cov=2, min_af=0.1, tile_size=93, tile_ov=45) {
    gr2df = function(gr) as_tibble(as.data.frame(gr)) %>%
        select(var_id, everything()) %>%
        arrange(order(gtools::mixedorder(var_id)))

    # changes peptide, is unique and is expressed
    subs = subset_context(res[! res$CONSEQUENCE %in% c("synonymous", "nonsense", "nostart")])
    alt_in_ref = function(a,r) grepl(as.character(a), as.character(r), fixed=TRUE)
    subs = subs[!mapply(alt_in_ref, a=subs$alt_prot, r=subs$ref_prot)]
    if ("rna_count" %in% colnames(S4Vectors::mcols(subs)) && !all(is.na(subs$rna_count)))
        subs = subs[!is.na(subs$rna_count) & subs$rna_count > 0] # & subs$rna_tpm > 0]
    subs = subs[subs$AF >= min_af & subs$cov_alt >= min_cov]
    subs = subs[!duplicated(subs$alt_prot)]

    # peptide is not contained within another peptide
    contained_in = function(i) any(grepl(ac2[i], ac2[-i], fixed=TRUE))
    ac = as.character(subs$alt_prot)
    ac2 = substr(ac, pmax(1, 1-subs$alt_shift/3), nchar(ac)-pmax(0, subs$alt_shift/3))
    any_con = sapply(seq_along(ac2), contained_in)
    subs = subs[!any_con]

    # name the variants
    mut_lab = ifelse(subs$CONSEQUENCE == "frameshift", "fs", subs$VARAA)
    pstarts = unname(sapply(subs$PROTEINLOC, function(p) p[[1]])) + subs$silent_start
    subs$mut_id = sprintf("%s_%s%i%s", subs$gene_name, subs$REFAA, pstarts, mut_lab)

    table(nchar(subs$alt_nuc))
    stopifnot(nchar(subs$alt_nuc) %% 3 == 0)

    # tile peptides to have max `tile_size` nt length
    tile_cDNA = function(p) {
        ntile = ceiling((nchar(p)-tile_ov) / (tile_size-tile_ov))
        starts = round(seq(0, nchar(p)-tile_size, length.out=ntile)/3) * 3 + 1
        stopifnot((starts[-length(starts)]+tile_size - starts[-1]) >= tile_ov)
        lapply(starts, function(s) substr(p, s, s+tile_size-1))
    }
    pep = with(subs, tibble(var_id=subs$var_id, mut_id=subs$mut_id,
                            gene_name=subs$gene_name, gene_id=subs$GENEID,
                            tx_id=subs$tx_name, ref=as.character(subs$ref_nuc),
                            alt=as.character(subs$alt_nuc))) %>%
        tidyr::pivot_longer(c(ref, alt), names_to="type", values_to="cDNA") %>%
        rowwise() %>% mutate(tiled = list(tile_cDNA(cDNA))) %>% ungroup() %>%
        mutate(n_tiles = sapply(tiled, length)) %>%
        tidyr::unnest(tiled) %>%
        mutate(pep_id = ifelse(type == "alt", mut_id, sub("([0-9]+)[a-zA-Z]+$", "\\1", mut_id)),
               tiled = unlist(tiled, use.names=FALSE),
               nt = nchar(tiled),
               peptide = as.character(translate(DNAStringSet(tiled), no.init.codon=TRUE))) %>%
        group_by(pep_id) %>%
            mutate(pep_id = if(n()>1) paste(pep_id, seq_along(pep_id), sep="-") else pep_id) %>%
        ungroup()

#    stopifnot(pep$type[duplicated(pep$tiled)] == "ref")
    stopifnot(!any(duplicated(pep$pep_id)))
    pep = pep[!duplicated(pep$tiled),]

    pep$BbsI_replaced = vcountPattern("GAAGAC", pep$tiled) + vcountPattern("GTCTTC", pep$tiled)
    pep$tiled = sapply(pep$tiled, remove_cutsite, site="GAAGAC", seed=178529, USE.NAMES=FALSE)
#    stopifnot(pep$peptide == as.character(translate(DNAStringSet(pep$tiled), no.init.codon=TRUE)))

    sv = list(
        `All Variants` = res %>% select(-CDSLOC) %>% gr2df() %>%
            dplyr::select(-tx_name, -(ref_nuc:alt_prot)) %>% distinct(),
        `Unique Protein-Coding` = gr2df(subs) %>% select(var_id, mut_id, everything()),
        `93 nt Peptides` = pep %>% select(var_id, mut_id, pep_id,
            gene_id:cDNA, n_tiles, BbsI_replaced, tiled, nt, peptide)
    )
    writexl::write_xlsx(sv, path=fname)
}
