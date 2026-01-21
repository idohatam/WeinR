#' Filter reads across all files in a LongReadQC object
#'
#' Applies per-read filters (mean Q-score, read length, and N-count) to each file
#' in a \code{LongReadQC} object and writes filtered reads to disk. A filtering
#' summary and output paths are stored in \code{qc_obj@metadata} for each file.
#'
#' @param qc_obj A \code{LongReadQC} object with per-read metrics already computed
#'   (e.g., via \code{QualMat()}).
#' @param OutDir Character(1). Output directory where filtered read files will be
#'   written. Created if it does not exist.
#' @param MinAvgQS Numeric(1). Minimum mean per-read Q-score required to keep a read.
#' @param MinLength Integer(1). Minimum read length required to keep a read.
#' @param MaxNumberNs Integer(1). Maximum number of \code{N} bases allowed per read.
#' @param OutFileType Character vector specifying output format(s). Supported values
#'   are \code{"fastq"}, \code{"fasta"}, and \code{"bam"}. Multiple formats may be
#'   requested.
#'
#' @return A \code{LongReadQC} object with per-file filtering results stored in
#'   \code{qc_obj@metadata[[file]]$filter_summary}.
#'
#' @seealso \code{\link{FilterLong}} for single-file filtering and
#'   \code{\link{QualMat}} for computing per-read metrics used by filtering.
#'
#' @examples
#' \dontrun{
#' qc_obj <- QualMat(qc_obj)  # ensure metrics exist first
#' qc_filt <- QualFilter(
#'   qc_obj,
#'   OutDir = "filtered_reads",
#'   MinAvgQS = 20,
#'   MinLength = 100,
#'   MaxNumberNs = 2,
#'   OutFileType = c("fastq", "fasta")
#' )
#' }
#'
#' @export
QualFilter <- function(qc_obj,
                       OutDir = ".",
                       MinAvgQS = 20,
                       MinLength = 100,
                       MaxNumberNs = 2,
                       OutFileType = c("fastq")) {
  
  message("Starting filtering workflow for ", length(qc_obj@files), " file(s)...")
  
  OutDir <- normalizePath(OutDir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  
  for (fpath in qc_obj@files) {
    
    message("\n-> Processing file: ", fpath)
    
    # Build single-file subobject
    single_qc <- qc_obj
    single_qc@files <- fpath
    single_qc@metrics[[fpath]] <- qc_obj@metrics[[fpath]]
    single_qc@summary_metrics <- qc_obj@summary_metrics[
      qc_obj@summary_metrics$file == fpath,
      , drop = FALSE
    ]
    
    # Run filtering
    filtered_qc <- FilterLong(
      qc_obj = single_qc,
      MinAvgQS = MinAvgQS,
      MinLength = MinLength,
      MaxNumberNs = MaxNumberNs,
      OutFileType = OutFileType,
      OutDir = OutDir
    )
    
    # If filtering produced no metadata, skip
    if (is.null(filtered_qc@metadata[[fpath]]$filter_summary)) {
      message("No filtering summary for: ", fpath)
      next
    }
    
    # Save ONLY this file's metadata
    qc_obj@metadata[[fpath]] <- filtered_qc@metadata[[fpath]]
  }
  
  message("\nFiltering complete for all files.")
  
  return(qc_obj)
}
