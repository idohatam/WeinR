#' Filter long-read sequences using QC metrics (internal)
#'
#' Internal helper that filters reads for a single input file referenced by a QC
#' object using per-read metrics (mean Q-score, read length, and N-count).
#' A per-file filter summary is stored in \code{qc_obj@metadata[[basename(file)]]$filter_summary}.
#' Optionally, an intermediate FASTQ is written for chaining and stored in
#' \code{attr(qc_obj, "._tmp_fastq")}.
#'
#' @param qc_obj A QC S4 object containing \code{@files} and \code{@metrics}.
#'   This function expects exactly one file in \code{qc_obj@files}.
#' @param MinAvgQS Numeric(1). Minimum mean per-read Q-score required to keep a read.
#' @param MinLength Integer(1). Minimum read length required to keep a read.
#' @param MaxNumberNs Integer(1). Maximum allowed number of \code{N} bases per read.
#' @param OutFileType Character vector specifying output format(s) to write when outputs
#'   are written (e.g., \code{"fastq"}). Supported values are \code{"fastq"},
#'   \code{"bam"}, and \code{"fasta"} (as implemented by \code{WriteReadOutputs()}).
#' @param OutDir Character(1). Output directory for filtered reads (created if needed).
#' @param WriteIntermediate Logical(1). If \code{TRUE}, write an intermediate FASTQ for chaining
#'   and store its path in \code{attr(qc_obj, "._tmp_fastq")}.
#' @param KeepIntermediates Logical(1). If \code{TRUE}, write filtered reads to \code{OutDir}
#'   (in the formats requested by \code{OutFileType}) and record written paths in the summary.
#' @param OutSuffix Character(1). Suffix appended to the output basename when writing filtered reads.
#'
#' @return The updated \code{qc_obj}. The filter summary is stored in
#'   \code{qc_obj@metadata[[basename(file)]]$filter_summary}. If \code{WriteIntermediate = TRUE},
#'   \code{attr(qc_obj, "._tmp_fastq")} contains the intermediate FASTQ path (or \code{NULL}).
#'
#' @details
#' If no reads pass filters, the function issues a warning, records a summary with
#' \code{reads_after = 0L} and \code{output_paths = character(0)}, and returns without writing outputs.
#'
#' @keywords internal

FilterLong <- function(qc_obj,
                       MinAvgQS = 20,
                       MinLength = 100,
                       MaxNumberNs = 2,
                       OutFileType = c("fastq"),
                       OutDir = ".",
                       WriteIntermediate = FALSE,
                       KeepIntermediates = FALSE,
                       OutSuffix = "filtered") {
  
  if (length(qc_obj@files) != 1L) {
    stop("FilterLong() expects a QC object with exactly one file in qc_obj@files.")
  }
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)
  
  # reset temp chaining attribute
  attr(qc_obj, "._tmp_fastq") <- NULL
  
  ## Extract file path and label
  fpath <- qc_obj@files[[1]]
  fname <- basename(fpath)
  
  ## Retrieve metrics for this file (keyed by basename)
  file_metrics <- qc_obj@metrics[[fname]]
  if (is.null(file_metrics)) {
    stop("Metrics not found for file: ", fpath)
  }
  
  ## Extract per-read metrics
  lengths  <- file_metrics$readLengths
  q_scores <- unlist(file_metrics$meanprQscore)
  n_counts <- unlist(file_metrics$Ncount)
  
  ## Load reads
  reads <- ImportFile(fpath)
  if (!inherits(reads, "QualityScaledDNAStringSet")) {
    stop("ImportFile() must return a QualityScaledDNAStringSet object.")
  }
  
  ## ---- consistency checks (invalid state -> error) ----
  n_reads <- length(reads)
  
  if (length(lengths) != n_reads ||
      length(q_scores) != n_reads ||
      length(n_counts) != n_reads) {
    stop(
      "FilterLong(): metric lengths do not match number of reads.\n",
      "  file: ", fname, "\n",
      "  reads: ", n_reads, "\n",
      "  readLengths: ", length(lengths), "\n",
      "  meanprQscore: ", length(q_scores), "\n",
      "  Ncount: ", length(n_counts)
    )
  }
  
  ## Compute keep vector (guaranteed aligned)
  keep <- (q_scores >= MinAvgQS) &
    (lengths  >= MinLength) &
    (n_counts <= MaxNumberNs)
  
  n_keep <- sum(keep)
  drop_frac <- (n_reads - n_keep) / n_reads
  
  ## ---- edge-case warnings ----
  if (n_keep == 0L) {
    warning("FilterLong(): no reads passed filters for: ", fname)
    
    qc_obj@metadata[[fname]]$filter_summary <- list(
      reads_before = n_reads,
      reads_after  = 0L,
      output_paths = character(0)
    )
    
    attr(qc_obj, "._tmp_fastq") <- NULL
    return(qc_obj)
  }
  
  if (drop_frac >= 0.50) {
    warning(
      "FilterLong(): filtering removed ",
      round(drop_frac * 100, 1), "% of reads for: ", fname,
      " (", n_keep, "/", n_reads, " kept)."
    )
  }
  
  filtered_reads <- reads[keep]
  
  out_paths <- character(0)
  if (isTRUE(KeepIntermediates)) {
    base_name <- paste0(RemoveExt(OriginalFPath), "_", OutSuffix)
    out_paths <- WriteReadOutputs(
      reads = final_reads,
      base_name = base_name,
      OutDir = OutDir,
      OutFileType = OutFileType
    )
  }
  qc_obj@metadata[[fname]]$filter_summary <- list(
    reads_before = n_reads,
    reads_after  = length(filtered_reads),
    output_paths = out_paths
  )
  
  ## Optional intermediate fastq for chaining
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

