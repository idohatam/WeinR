#' Filter long-read sequences using QC metrics (internal)
#'
#' Internal helper that filters reads for a single file represented inside a QC
#' object using per-read metrics (mean Q-score, read length, and N-count).
#' Reads that pass the filters are written to disk, and a filter summary is stored
#' in \code{qc_obj@metadata[[basename(file)]]}.
#'
#' @param qc_obj A QC S4 object containing \code{@files} and \code{@metrics}.
#'   This function expects exactly one file in \code{qc_obj@files}.
#' @param MinAvgQS Numeric(1). Minimum mean per-read Q-score required to keep a read.
#' @param MinLength Integer(1). Minimum read length required to keep a read.
#' @param MaxNumberNs Integer(1). Maximum allowed number of \code{N} bases per read.
#' @param OutFileType Character vector specifying output format(s).
#'   Supported values are \code{"fastq"}, \code{"bam"}, and \code{"fasta"}.
#'   Multiple formats may be requested.
#' @param OutDir Character(1). Output directory for filtered reads.
#'   Created if it does not exist.
#' @param WriteIntermediate Logical(1). If \code{TRUE}, an intermediate FASTQ file
#'   is written and its path is stored in
#'   \code{attr(qc_obj, "._tmp_fastq")} for chaining steps in downstream workflows.
#'
#' @return The updated \code{qc_obj}. If no reads pass the filters, the input
#'   object is returned unchanged and no output files are written.
#'
#' @details
#' Per-read QC metrics are retrieved from \code{qc_obj@metrics} using the
#' basename of the input file as the key. Filtered reads are written using
#' \code{WriteReadOutputs()}, and a summary of filtering results is stored in
#' \code{qc_obj@metadata[[basename(file)]]}.
#'
#' When \code{WriteIntermediate = TRUE}, a temporary FASTQ file containing the
#' filtered reads is written to a session-specific directory under
#' \code{tempdir()}, and its path is attached to the QC object for use by
#' subsequent processing steps.
#'
#' @examples
#' \dontrun{
#' qc_obj2 <- FilterLong(
#'   qc_obj,
#'   MinAvgQS = 20,
#'   MinLength = 100,
#'   MaxNumberNs = 2,
#'   OutFileType = c("fastq", "fasta"),
#'   OutDir = tempdir(),
#'   WriteIntermediate = TRUE
#' )
#' }
#'
#' @keywords internal
FilterLong <- function(qc_obj,
                       MinAvgQS = 20,
                       MinLength = 100,
                       MaxNumberNs = 2,
                       OutFileType = c("fastq"),
                       OutDir = ".",
                       WriteIntermediate = FALSE) {
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  
  attr(qc_obj, "._tmp_fastq") <- NULL
  
  ## Extract file path and label
  fpath <- qc_obj@files[[1]]
  fname <- basename(fpath)
  
  ## Retrieve metrics for this file
  file_metrics <- qc_obj@metrics[[fname]]
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
  
  
  qc_obj@metadata[[fname]]$filter_summary <- list(
    reads_before = length(reads),
    reads_after  = length(filtered_reads),
    output_paths = out_paths
  )
  
  tmp_fastq <- NULL
  if (isTRUE(WriteIntermediate) && length(filtered_reads) > 0L) {
    
    stem <- RemoveExt(basename(fpath))
    tmp_base <- paste0(stem, "_filtered_tmp_", Sys.getpid(), "_", sample.int(1e9, 1))
    
    tmp_dir <- file.path(tempdir(), "longR_intermediate")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    
    tmp_paths <- WriteReadOutputs(
      reads = filtered_reads,
      base_name = tmp_base,
      OutDir = tmp_dir,
      OutFileType = "fastq",
      Verbose = FALSE
    )
    
    tmp_fastq <- tmp_paths[[1]]
  }
  

  
  attr(qc_obj, "._tmp_fastq") <- tmp_fastq
  
  base::return(qc_obj)
}

