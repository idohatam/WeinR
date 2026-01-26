#' Import a single sequence file as a QualityScaledDNAStringSet
#'
#' Internal helper that reads one input file (FASTQ/FASTQ.GZ/BAM) and returns a
#' \code{Biostrings::QualityScaledDNAStringSet}. FASTQ inputs are pre-checked with
#' \code{FastqCheck()} before reading. BAM inputs are opened with
#' \code{Rsamtools::BamFile()} and read via \code{Rsamtools::scanBam()} then
#' converted to a quality-scaled set.
#'
#' @param filePath Character scalar. Path to a single input file.
#' @param MinNumReads Integer(1). Minimum number of reads required to keep the file.
#'   Must be a single positive integer (>= 1). If fewer reads are present, the
#'   function warns and returns \code{NULL}.
#'
#' @return A \code{Biostrings::QualityScaledDNAStringSet} on success, or \code{NULL}
#'   if the number of reads is less than \code{MinNumReads}.
#'
#' @details
#' File type is inferred by \code{CheckFile()}. For FASTQ/FASTQ.GZ inputs, the file
#' is first validated with \code{FastqCheck(filePath, n_records = 100L)}; failures
#' are rethrown as \dQuote{FASTQ check failed: ...}. Reading failures from
#' \code{Biostrings::readQualityScaledDNAStringSet()} are rethrown as
#' \dQuote{FASTQ read failed (corrupted/malformed): ...}.
#'
#' For BAM inputs, the file is opened and read with \code{Rsamtools}. Errors while
#' opening or scanning are rethrown as \dQuote{BAM open failed ...} or
#' \dQuote{BAM read failed ...}. The function also validates required BAM fields
#' (\code{seq} and \code{qual}) and checks that their lengths match. If available
#' and aligned in length, BAM \code{qname} values are used as read names.
#'
#' This function is intended for internal use and is not part of the public API.
#'
#' @examples
#' \dontrun{
#' # FASTQ example: create a tiny, valid FASTQ file (1 read)
#' fq <- tempfile(fileext = ".fastq")
#' writeLines(c("@r1", "ACGT", "+", "IIII"), fq)
#' x <- WeinR:::ImportFile(fq)
#' length(x)
#'
#' # MinNumReads example: warn + return NULL if threshold not met
#' y <- WeinR:::ImportFile(fq, MinNumReads = 2L)
#' is.null(y)
#'
#' # MinNumReads validation example: errors on invalid values
#' WeinR:::ImportFile(fq, MinNumReads = 0)
#' }
#'
#' @keywords internal
ImportFile <- function(filePath,
                       MinNumReads = 1L) {
  
  
  # validate MinNumReads
  if (is.null(MinNumReads) || !is.numeric(MinNumReads) ||
      length(MinNumReads) != 1L || is.na(MinNumReads) ||
      MinNumReads < 1 || MinNumReads %% 1 != 0) {
    stop("MinNumReads must be a single positive integer (>= 1).")
  }
  MinNumReads <- as.integer(MinNumReads)
  
  infile <- CheckFile(filePath)
  
  Output <- if (infile %in% c("fastq", "fastq.gz")) {
    
    tryCatch(FastqCheck(filePath, n_records = 100L),
             error = function(e) stop("FASTQ check failed: ", conditionMessage(e)))
    
    tryCatch(
      Biostrings::readQualityScaledDNAStringSet(filePath),
      error = function(e) {
        stop(sprintf("FASTQ read failed (corrupted/malformed): %s", conditionMessage(e)))
      }
    )
    
  } else if (infile == "bam") {
    
    bf <- Rsamtools::BamFile(filePath)
    on.exit(try(Rsamtools::close(bf), silent = TRUE), add = TRUE)
    
    tryCatch(Rsamtools::open(bf),
             error = function(e) {
               stop(sprintf("BAM open failed (corrupted/invalid): %s", conditionMessage(e)))
             })
    
    bam <- tryCatch(
      Rsamtools::scanBam(bf)[[1]],
      error = function(e) {
        stop(sprintf("BAM read failed (corrupted/invalid): %s", conditionMessage(e)))
      }
    )
    
    if (is.null(bam$seq) || is.null(bam$qual)) {
      stop("BAM missing required fields 'seq' and/or 'qual'.")
    }
    if (length(bam$seq) != length(bam$qual)) {
      stop("BAM 'seq' and 'qual' lengths do not match.")
    }
    
    out <- Biostrings::QualityScaledDNAStringSet(bam$seq, bam$qual)
    
    if (length(out) > 0L && !is.null(bam$qname) && length(bam$qname) == length(out)) {
      names(out) <- bam$qname
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
  
  Output
}
