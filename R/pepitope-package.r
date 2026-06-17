#' Pepitope: peptide epitopes from reference genome and variant (VCF) file
#'
#' Given a reference genome and VCF file, this package will provide
#' the upstream/downstream peptide context of a variant. It will generate a
#' summary report for protein-coding variants including the reference and
#' mutated allele, read coverage, amino acid sequence, and other information.
#' It can also be used to remove restriction sites from cDNA, alongside other
#' helper functions.
#'
#' @return Package documentation for pepitope
#' @name pepitope
#' @docType package
#' @import dplyr
#' @import ggplot2
#' @importFrom Rcpp evalCpp
#' @useDynLib pepitope, .registration=TRUE
"_PACKAGE"
