% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_fusions.r
\name{filter_fusions}
\alias{filter_fusions}
\title{Filter a fusion VRanges object by number of reads and tools}
\usage{
filter_fusions(vr, min_reads = NULL, min_split_reads = NULL, min_tools = NULL)
}
\arguments{
\item{vr}{A VRanges object with RNA fusions from readVcfAsRanges}

\item{min_reads}{The minimum number of linked read support for a fusion}

\item{min_split_reads}{The minimum number of split read support for a fusion}

\item{min_tools}{The minimum number of tools that identify a fusion}
}
\value{
A filtered VRanges object
}
\description{
Filter a fusion VRanges object by number of reads and tools
}
