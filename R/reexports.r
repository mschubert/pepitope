#' Objects exported from other packages
#'
#' These objects are imported from other packages. Follow the links
#' below to see their documentation.
#'
#' \describe{
#'   \item{GenomeInfoDb}{\code{\link[Seqinfo:seqinfo]{genome}},
#'     \code{\link[Seqinfo:seqinfo]{seqinfo}},
#'     \code{\link[Seqinfo:seqinfo]{seqlevels}},
#'     \code{\link[GenomeInfoDb:seqlevelsStyle]{seqlevelsStyle()}},
#'     \code{\link[GenomeInfoDb:seqlevelsStyle<-]{seqlevelsStyle<-()}}}
#'   \item{SummarizedExperiment}{\code{\link[SummarizedExperiment:assay]{assay()}},
#'     \code{\link[SummarizedExperiment:colData]{colData()}},
#'     \code{\link[SummarizedExperiment:rowData]{rowData()}}}
#'   \item{VariantAnnotation}{\code{\link[VariantAnnotation:VRanges-class]{readVcfAsVRanges()}}}
#' }
#'
#' @docType import
#' @name reexports
#' @keywords internal
NULL

#' @rdname reexports
#' @importFrom SummarizedExperiment assay
#' @export
SummarizedExperiment::assay

#' @rdname reexports
#' @importFrom SummarizedExperiment colData
#' @export
SummarizedExperiment::colData

#' @rdname reexports
#' @importFrom SummarizedExperiment rowData
#' @export
SummarizedExperiment::rowData

#' @rdname reexports
#' @importFrom GenomeInfoDb genome
#' @export
GenomeInfoDb::genome

#' @rdname reexports
#' @importFrom GenomeInfoDb seqinfo
#' @export
GenomeInfoDb::seqinfo

#' @rdname reexports
#' @importFrom GenomeInfoDb seqlevels
#' @export
GenomeInfoDb::seqlevels

#' @rdname reexports
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @export
GenomeInfoDb::seqlevelsStyle

#' @rdname reexports
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @export
GenomeInfoDb::`seqlevelsStyle<-`

#' @rdname reexports
#' @importFrom VariantAnnotation readVcfAsVRanges
#' @export
VariantAnnotation::readVcfAsVRanges
