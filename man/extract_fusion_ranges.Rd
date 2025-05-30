% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_fusions.r
\name{extract_fusion_ranges}
\alias{extract_fusion_ranges}
\title{Convert a fusion VRanges object to left (5') and right (3') GRanges objects}
\usage{
extract_fusion_ranges(vr)
}
\arguments{
\item{vr}{A VRanges object with RNA fusions from `readVcfAsVRanges`}
}
\value{
A list of the 5' and 3' GRanges objects
}
\description{
Convert a fusion VRanges object to left (5') and right (3') GRanges objects
}
\keyword{internal}
