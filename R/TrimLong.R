#' Trim reads from the start and/or end (internal)
#'
#' Internal helper that trims reads for a single input file using \code{Start} and
#' \code{End} offsets and writes trimmed reads to disk. A trim summary is stored in
#' \code{qc_obj@metadata[[basename(original_file)]]}.
#'
#' @param qc_obj A QC S4 object containing \code{@files}. This function expects
#'   exactly one file in \code{qc_obj@files} (the "original" file key for metadata).
#' @param Start Integer(1) or \code{NULL}. Number of bases to remove from the start
#'   (5') of each read. If \code{NULL}, no start trimming is applied.
#' @param End Integer(1) or \code{NULL}. Number of bases to remove from the end
#'   (3') of each read. If \code{NULL}, no end trimming is applied.
#' @param FilePath Character(1) or \code{NULL}. Optional file path to trim. If
#'   \code{NULL}, uses the first file in \code{qc_obj@files}. This is used to allow
#'   trimming of intermediate files while still recording results under the original
#'   input file.
#' @param OutDir Character(1). Output directory. Created if it does not exist.
#' @param OutFile Character(1) or \code{NULL}. Optional output filename stem. If
#'   provided, \code{RemoveExt(OutFile)} is used as the base name; otherwise the
#'   original input file stem is used.
#' @param OutFileType Character vector specifying output format(s). Supported values
#'   are \code{"fastq"}, \code{"fasta"}, and \code{"bam"} (as implemented by
#'   \code{WriteReadOutputs()}). Multiple formats may be requested.
#'
#' @return The updated \code{qc_obj}. A trim summary is stored in
#'   \code{qc_obj@metadata[[basename(original_file)]]$trim_summary}.
#'
#' @details
#' Trimming is performed using \code{IRanges::narrow()}. Reads that would become
#' invalid (end < start) are dropped. If no reads can be trimmed under the provided
#' \code{Start}/\code{End} values, the input \code{qc_obj} is returned unchanged.
#'
#' @examples
#' \dontrun{
#' qc_obj2 <- TrimLong(qc_obj, Start = 10, End = 5, OutDir = tempdir())
#' }
#'
#' @keywords internal
TrimLong <- function(qc_obj,
                     Start = NULL,
                     End = NULL,
                     FilePath = NULL,
                     OutDir = ".",
                     OutFile = NULL,
                     OutFileType = c("fastq")) {
  
  # Choose which file to trim
  fpath <- if (!is.null(FilePath)) FilePath else qc_obj@files[[1]]
  OriginalFPath <- qc_obj@files[[1]]
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  
  # Import reads
  reads <- ImportFile(fpath)
  
  # Compute trimming
  n <- Biostrings::width(reads)
  trim_start <- rep(if (!is.null(Start)) Start + 1 else 1, length(n))
  trim_end   <- if (!is.null(End)) n - End else n
  
  keep <- trim_end >= trim_start
  
  if (!any(keep)) {
    message("TrimLong(): no reads can be trimmed with these Start/End values.")
    return(qc_obj)
  }
  
  trimmed_reads <- IRanges::narrow(
    reads[keep],
    start = trim_start[keep],
    end   = trim_end[keep]
  )
  #trimmed_reads <- IRanges::narrow(reads, start = trim_start, end = trim_end)
  
  # Output naming
  base_name <- if (!is.null(OutFile)) RemoveExt(OutFile) else RemoveExt(OriginalFPath)
  
  out_paths <- WriteReadOutputs(
    reads = trimmed_reads,
    base_name = paste0(base_name, "_trimmed"),
    OutDir = OutDir,
    OutFileType = OutFileType
  )
  
  # metadata storage (per file path)
  qc_obj@metadata[[basename(OriginalFPath)]]$trim_summary <- list(
    reads_before = length(reads),
    reads_after = length(trimmed_reads),
    output_paths = out_paths
  )
  
  base::return(qc_obj)
}

