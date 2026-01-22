#' Import a single sequence file as a QualityScaledDNAStringSet
#'
#' Internal helper that reads one input file (FASTQ/FASTQ.GZ/BAM) and returns a
#' `Biostrings::QualityScaledDNAStringSet`. BAM files are read with
#' `Rsamtools::scanBam()` and converted to a quality-scaled set.
#'
#' @param filePath Character scalar. Path to a single input file.
#' @param MinNumReads Integer(1). Minimum number of reads required to keep the file.
#'   If fewer reads are present, the function warns and returns `NULL`.
#'
#' @return A `Biostrings::QualityScaledDNAStringSet` on success, or `NULL` if the
#'   number of reads is less than `MinNumReads`.
#'
#' @details
#' File type is inferred by `CheckFile()`. Reading errors are rethrown with a more
#' informative message (e.g., malformed FASTQ or invalid BAM).
#'
#' This function is intended for internal use and is not part of the public API.
#'
#' @examples
#' \dontrun{
#' # FASTQ example: create a tiny, valid FASTQ file
#' fq <- tempfile(fileext = ".fastq")
#' writeLines(c("@r1", "ACGT", "+", "IIII"), fq)
#' x <- WeinR:::ImportFile(fq)
#' length(x)
#'
#' # MinNumReads example: warn + return NULL if threshold not met
#' y <- WeinR:::ImportFile(fq, MinNumReads = 2L)
#' is.null(y)
#' }
#'
#' @keywords internal
ImportFile <- function(filePath, MinNumReads = 1L) {
  
  # Validate MinNumReads (empty files will be skipped since MinNumReads >= 1)
  if (is.null(MinNumReads) || !is.numeric(MinNumReads) ||
      length(MinNumReads) != 1L || is.na(MinNumReads) ||
      MinNumReads < 1 || MinNumReads %% 1 != 0) {
    stop("MinNumReads must be a single positive integer (>= 1).")
  }
  MinNumReads <- as.integer(MinNumReads)
  
  infile <- CheckFile(filePath)
  
  # Read file (catch corruption / unreadable file errors with clearer messaging)
  Output <- if (infile %in% c("fastq", "fastq.gz")) {
    
    tryCatch(
      Biostrings::readQualityScaledDNAStringSet(filePath),
      error = function(e) {
        stop(sprintf("FASTQ read failed (corrupted/malformed): %s", conditionMessage(e)))
      }
    )
    
  } else if (infile == "bam") {
    
    bam <- tryCatch(
      Rsamtools::scanBam(filePath)[[1]],
      error = function(e) {
        stop(sprintf("BAM read failed (corrupted/invalid): %s", conditionMessage(e)))
      }
    )
    
    out <- Biostrings::QualityScaledDNAStringSet(bam$seq, bam$qual)
    
    # Attach read names if available (only warn if there are reads)
    if (length(out) > 0L) {
      if (!is.null(bam$qname) && length(bam$qname) == length(out)) {
        names(out) <- bam$qname
      } else {
        warning("BAM qname missing or length mismatch; read names not attached.")
      }
    }
    
    out
  }
  
  # Minimum reads threshold (edge case: warn + skip)
  n_reads <- length(Output)
  if (n_reads < MinNumReads) {
    warning(sprintf(
      "File has %d reads (< MinNumReads = %d). File skipped: %s",
      n_reads, MinNumReads, filePath
    ))
    return(NULL)
  }
  
 return (Output)
}
