#' Trim reads from the start and/or end (internal)
#'
#' Internal helper that trims reads for a single input file using \code{Start} and/or
#' \code{End} offsets and writes trimmed reads to disk. A trim summary is stored in
#' \code{qc_obj@metadata[[basename(original_file)]]$trim_summary}.
#'
#' Edge-case behavior:
#' \itemize{
#'   \item Warns if trimming removes \code{>= 50\%} of reads.
#'   \item Warns if \code{0} reads remain after trimming; no files are written.
#' }
#'
#' @param qc_obj A QC S4 object containing \code{@files}. This function expects
#'   exactly one file in \code{qc_obj@files} (the "original" file key for metadata).
#' @param Start Integer(1) or \code{NULL}. Number of bases to remove from the start
#'   (5') of each read. If \code{NULL}, no start trimming is applied.
#' @param End Integer(1) or \code{NULL}. Number of bases to remove from the end
#'   (3') of each read. If \code{NULL}, no end trimming is applied.
#' @param FilePath Character(1) or \code{NULL}. Optional file path to trim. If
#'   \code{NULL}, uses the first file in \code{qc_obj@files}. This allows trimming of
#'   intermediate files while still recording results under the original input file.
#' @param OutDir Character(1). Output directory. Created if it does not exist.
#' @param OutFileType Character vector specifying output format(s). Supported values
#'   are \code{"fastq"}, \code{"fasta"}, and \code{"bam"} (as implemented by
#'   \code{WriteReadOutputs()}). Multiple formats may be requested.
#' @param OutSuffix Character(1). Suffix appended to the output basename when writing
#'   trimmed reads (e.g., \code{"trimmed"}).
#'
#' @return The updated \code{qc_obj}. A trim summary is stored in
#'   \code{qc_obj@metadata[[basename(original_file)]]$trim_summary}.
#'
#' @details
#' Trimming is performed using \code{IRanges::narrow()}. Reads that would become
#' invalid (\code{end < start}) are dropped before trimming. If \code{0} reads remain,
#' a warning is issued, the summary is recorded with \code{reads_after = 0L}, and no
#' output files are written.
#'
#' @keywords internal
TrimLong <- function(qc_obj,
                     Start = NULL,
                     End = NULL,
                     FilePath = NULL,
                     OutDir = ".",
                     OutFileType = c("fastq"),
                     OutSuffix = "trimmed") {
  
  if (length(qc_obj@files) != 1L) {
    stop("TrimLong() expects a QC object with exactly one file in qc_obj@files.")
  }
  
  # Validate Start/End
  .chk_int1 <- function(x, nm) {
    if (is.null(x)) return(invisible(TRUE))
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x %% 1 != 0) {
      stop("TrimLong(): ", nm, " must be NULL or a single non-negative integer.")
    }
    invisible(TRUE)
  }
  .chk_int1(Start, "Start")
  .chk_int1(End, "End")
  
  # Choose which file to trim
  fpath <- if (!is.null(FilePath)) FilePath else qc_obj@files[[1]]
  OriginalFPath <- qc_obj@files[[1]]
  key <- basename(OriginalFPath)  # metadata key always uses original input
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)
  
  # Import reads
  reads <- ImportFile(fpath)
  if (!inherits(reads, "QualityScaledDNAStringSet")) {
    stop("ImportFile() must return a QualityScaledDNAStringSet object.")
  }
  
  n_reads <- length(reads)
  if (n_reads == 0L) {
    warning("TrimLong(): input contains 0 reads for: ", key)
    qc_obj@metadata[[key]]$trim_summary <- list(
      reads_before = 0L,
      reads_after  = 0L,
      output_paths = character(0)
    )
    return(qc_obj)
  }
  
  # Compute trimming bounds (per-read lengths)
  n <- Biostrings::width(reads)
  Start_i <- if (!is.null(Start)) as.integer(Start) else NULL
  End_i   <- if (!is.null(End))   as.integer(End)   else NULL
  
  trim_start <- rep(if (!is.null(Start_i)) Start_i + 1L else 1L, n_reads)
  trim_end   <- if (!is.null(End_i)) n - End_i else n
  
  keep <- trim_end >= trim_start
  n_after <- sum(keep)
  drop_frac <- (n_reads - n_after) / n_reads
  
  # Edge-case: 0 reads remain -> warn, store summary, write nothing
  if (n_after == 0L) {
    warning("TrimLong(): no reads remain after trimming for: ", key)
    
    qc_obj@metadata[[key]]$trim_summary <- list(
      reads_before = n_reads,
      reads_after  = 0L,
      output_paths = character(0)
    )
    
    return(qc_obj)
  }
  
  # Edge-case: large drop -> warn
  if (drop_frac >= 0.50) {
    warning(
      "TrimLong(): trimming removed ",
      round(drop_frac * 100, 1), "% of reads for: ", key,
      " (", n_after, "/", n_reads, " kept)."
    )
  }
  
  # Apply trimming
  trimmed_reads <- IRanges::narrow(
    reads[keep],
    start = trim_start[keep],
    end   = trim_end[keep]
  )
  
 
  # Output naming
  base_name <- paste0(RemoveExt(OriginalFPath))
  
  out_paths <- WriteReadOutputs(
    reads = trimmed_reads,
    base_name = paste0(base_name, "_", OutSuffix),
    OutDir = OutDir,
    OutFileType = OutFileType
  )
  
  
  # Store summary under original key
  qc_obj@metadata[[key]]$trim_summary <- list(
    reads_before = n_reads,
    reads_after  = length(trimmed_reads),
    output_paths = out_paths
  )
  
  base::return(qc_obj)
}

