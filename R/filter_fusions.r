#' Filter a fusion VRanges object by number of reads and tools
#'
#' @param vr    A VRanges object with RNA fusions from readVcfAsRanges
#' @param min_reads  The minimum number of linked read support for a fusion
#' @param min_split_reads  The minimum number of split read support for a fusion
#' @param min_tools  The minimum number of tools that identify a fusion
#' @return      A filtered VRanges object
#'
#' @export
filter_fusions = function(vr, min_reads=NULL, min_split_reads=NULL, min_tools=NULL) {
    # GT: ref allele/i alt allele (always './1'?)
    # DV: split reads
    # RV: discordant mates supporting translocation
    # FFPM: fusion fragments per million RNA fragments
    if (!is.null(min_split_reads))
        vr = vr[vr$DV >= min_split_reads]
    if (!is.null(min_reads))
        vr = vr[vr$DV + vr$RV >= min_reads]
    if (!is.null(min_tools))
        vr = vr[unlist(vr$TOOL_HITS) >= min_tools]

    vr
}
