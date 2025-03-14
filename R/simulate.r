#' Simulate mutations and write them into a VCF file
#'
#' @param fname File name of the VCF to save
#' @importFrom dplyr transmute
#' @keywords internal
sim_vcf = function(fname) {
    # start from the "pan-can hotspots", then add some specific test cases
    ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
    asm = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38

    hotspots = readr::read_tsv("hotspots.txt") |>
        transmute(gene_name = Gene,
                  residue = sub("[0-9]+$", "", Residue),
                  position = `Amino Acid Position`,
                  variants = Variants)

    # mapping gene names, ids, and transcripts
    genes = genes(ens106) |>
        as.data.frame() |> as_tibble() |>
        dplyr::filter(gene_biotype == "protein_coding",
                      gene_name %in% hotspots$gene_name) |>
        dplyr::select(seqnames, gene_id, gene_name) |>
        distinct()

    # get transcript+protein annotations for hotspot gene names
    hstx = transcripts(ens106) |>
        as.data.frame() |> as_tibble() |>
        dplyr::filter(tx_biotype == "protein_coding",
                      seqnames %in% c(1:22,'X')) |>
        dplyr::select(seqnames, gene_id, tx_id, tx_is_canonical) |>
        inner_join(genes) |>
        inner_join(as.data.frame(proteins(ens106))) |>
        inner_join(hotspots, relationship="many-to-many") |>
        dplyr::filter(residue == substr(protein_sequence, position, position)) |>
        group_by(gene_name, position) |>
            slice_max(tx_is_canonical, n=1, with_ties=FALSE) |>
        ungroup()
    #stopifnot(setequal(hotspots$gene_name, hstx$gene_name)) # PIK3R1 missing

    # convert protein to genomic positions
    locs = IRanges::IRanges(start=as.integer(hstx$position), end=as.integer(hstx$position))
    names(locs) = hstx$protein_id
    res = do.call(c, unname(proteinToGenome(locs, ens106)))
    res$ref_nuc = BSgenome::getSeq(asm, res) # strand="+" does not work and not error
    stopifnot(Biostrings::translate(res$ref_nuc, no.init.codon=TRUE) == hstx$residue)
    stopifnot(hstx$tx_id == res$tx_id)
    res$ref_nuc[strand(res) == "-"] = Biostrings::reverseComplement(res$ref_nuc[strand(res) == "-"])
    res$tx_id = NULL

    # convert possible variants to long format
    vars = as.data.frame(res) |>
        cbind(hstx |> dplyr::select(tx_id, gene_id, gene_name, residue, position, variants)) |>
        mutate(ref_nuc = as.character(res$ref_nuc),
               var = strsplit(variants, "\\|")) |>
        tidyr::unnest(var) |>
        mutate(n_pts = as.integer(sub("[A-Z*]+:", "", var)),
               var = sub(":[0-9]+", "", var)) |>
        dplyr::filter(residue != var) |>
        dplyr::select(-variants)

    # use any codon that fits alternative variants
    codons = Biostrings::getGeneticCode()
    codons = codons[!duplicated(codons)] # this always takes first, use minimal change?
    codons = setNames(names(codons), codons)
    alt_nuc = Biostrings::DNAStringSet(codons[vars$var])
    alt_nuc[vars$strand == "-"] = Biostrings::reverseComplement(alt_nuc[vars$strand == "-"])

    # convert to VRanges that pepitope can handle
    vr = VariantAnnotation::VRanges(sampleNames = "tumor",
                                    seqnames = vars$seqnames,
                                    ranges = IRanges(vars$start, vars$end),
                                    ref = vars$ref_nuc,
                                    alt = alt_nuc,
                                    n_pts = vars$n_pts)

    VariantAnnotation::writeVcf(vr, fname)
}

#' Simulate sequencing data and write them to a FASTQ file
#'
#' @keywords internal
sim_seq = function() {
    # multiple samples: noise + 2 clear drop-outs

    # merge with barcodes

    # use biostrings to write fastq
}
