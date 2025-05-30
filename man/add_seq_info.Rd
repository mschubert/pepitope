% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_fusions.r
\name{add_seq_info}
\alias{add_seq_info}
\title{Adds sequence information to break transcripts}
\usage{
add_seq_info(gr, cds_break, asm, txdb, tx)
}
\arguments{
\item{gr}{GenomicRanges object of break location}

\item{cds_break}{A list of transcripts overlapping break from `cds_by_break`}

\item{asm}{A Genome sequence package object, eg. ::BSgenome.Hsapiens.NCBI.GRCh38}

\item{txdb}{A transcription database, eg. AnnotationHub()[["AH100643"]]}

\item{tx}{A list of transcripts obtained from `transcripts(txdb)`}
}
\value{
A DataFrame with sequence information
}
\description{
Adds sequence information to break transcripts
}
\keyword{internal}
