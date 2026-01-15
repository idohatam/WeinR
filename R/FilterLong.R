#' Filter long-read sequences using QC metrics (internal)
#'
#' Internal helper that filters reads for a single file represented inside a QC
#' object using per-read metrics (mean Q-score, read length, and N-count). Reads
#' that pass are written to disk, and a filter summary is stored in
#' \code{qc_obj@metadata[[file]]}.
#'
#' @param qc_obj A QC S4 object containing \code{@files} and \code{@metrics}.
#' @param MinAvgQS Numeric(1). Minimum mean per-read Q-score to keep a read.
#' @param MinLength Integer(1). Minimum read length to keep a read.
#' @param MaxNumberNs Integer(1). Maximum allowed number of N bases per read.
#' @param OutFileType Character vector specifying output format(s).
#'   Supported values are \code{"fastq"}, \code{"bam"}, and \code{"fasta"}.
#'   Multiple formats may be requested.
#' @param OutDir Character(1). Output directory. Created if it does not exist.
#'
#' @return The updated \code{qc_obj}. If no reads pass filters, the input
#'   \code{qc_obj} is returned unchanged.
#'
#' @details
#' Output file types are validated using \code{match.arg()}. Files are written
#' using \code{WriteReadOutputs()}, and a summary of filtering results is stored
#' in \code{qc_obj@metadata}.
#'
#' @examples
#' \dontrun{
#' qc_obj2 <- FilterLong(
#'   qc_obj,
#'   MinAvgQS = 20,
#'   MinLength = 100,
#'   MaxNumberNs = 2,
#'   OutFileType = c("fastq", "fasta"),
#'   OutDir = tempdir()
#' )
#' }
#'
#' @keywords internal
FilterLong <- function(qc_obj,
                       MinAvgQS = 20,
                       MinLength = 100,
                       MaxNumberNs = 2,
                       OutFileType = c("fastq"),
                       OutDir = ".") {
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  
  ## Extract file path and label
  fpath <- qc_obj@files[[1]]
  fname <- basename(fpath)
  
  ## Retrieve metrics for this file
  file_metrics <- qc_obj@metrics[[fpath]]
  if (is.null(file_metrics))
    stop("Metrics not found for file: ", fpath)
  
  ## Extract per-read metrics
  lengths  <- file_metrics$readLengths
  q_scores <- unlist(file_metrics$meanprQscore)
  n_counts <- unlist(file_metrics$Ncount)
  
  keep <- (q_scores >= MinAvgQS) &
    (lengths  >= MinLength) &
    (n_counts <= MaxNumberNs)
  
  ## Load reads
  reads <- ImportFile(fpath)
  
  if (!any(keep)) {
    message("No reads passed filters for: ", fname)
    return(qc_obj)
  }
  
  filtered_reads <- reads[keep]
  
  ## Write output
  base_name <- paste0(RemoveExt(fpath), "_filtered")
  
  out_paths <- WriteReadOutputs(
    reads = filtered_reads,
    base_name = base_name,
    OutDir = OutDir,
    OutFileType = OutFileType
  )
  
  
  qc_obj@metadata[[fpath]] <- list(
    filter_summary = list(
      reads_before = length(reads),
      reads_after  = length(filtered_reads),
      output_paths = out_paths,
      filter_date = Sys.time()
    )
  )
  
  base::return(qc_obj)
}

