#' Sniff a FASTQ/FASTQ.GZ file for basic structural validity
#'
#' Internal helper that reads the first \code{n_records} FASTQ records and checks
#' basic file structure:
#' \itemize{
#'   \item Records occur in groups of four lines
#'   \item Header lines (line 1 of each record) start with \code{"@"}
#'   \item Plus lines (line 3 of each record) start with \code{"+"}
#' }
#'
#' This function is intended as a lightweight corruption check and does not
#' validate sequence/quality lengths or encoding.
#'
#' @param filePath Character(1). Path to a \code{.fastq}, \code{.fq}, or gzipped
#'   equivalent (\code{.gz}).
#' @param n_records Integer(1). Number of FASTQ records to inspect. Must be a
#'   single positive integer (>= 1).
#'
#' @return Invisibly returns \code{TRUE} on success. Errors if the file appears
#'   malformed. An empty file is treated as valid and returns \code{TRUE}.
#'
#' @details
#' The function reads at most \code{n_records * 4} lines from the file. If the
#' number of lines read is not a multiple of four, or if header or plus lines do
#' not conform to FASTQ conventions, an error is raised with an informative
#' message.
#'
#' @keywords internal
FastqCheck <- function(filePath, n_records = 100L) {
  
  if (is.null(n_records) || !is.numeric(n_records) || length(n_records) != 1L ||
      is.na(n_records) || n_records < 1 || n_records %% 1 != 0) {
    stop("n_records must be a single positive integer (>= 1).")
  }
  n_records <- as.integer(n_records)
  
  n_lines <- n_records * 4L
  
  con <- NULL
  on.exit({
    if (!is.null(con)) try(close(con), silent = TRUE)
  }, add = TRUE)
  
  con <- if (grepl("\\.gz$", filePath, ignore.case = TRUE)) {
    gzfile(filePath, open = "rt")
  } else {
    file(filePath, open = "rt")
  }
  
  x <- readLines(con, n = n_lines, warn = FALSE)
  
  # Empty file is not corruption here
  if (length(x) == 0L) return(invisible(TRUE))
  
  # Must end on record boundary
  if ((length(x) %% 4L) != 0L) {
    stop("FASTQ appears malformed: incomplete record (lines not multiple of 4).")
  }
  
  # Header and plus-line checks
  h <- x[seq(1L, length(x), by = 4L)]
  p <- x[seq(3L, length(x), by = 4L)]
  
  if (any(!startsWith(h, "@"))) {
    stop("FASTQ appears malformed: header line(s) not starting with '@'.")
  }
  if (any(!startsWith(p, "+"))) {
    stop("FASTQ appears malformed: plus line(s) not starting with '+'.")
  }
  
  invisible(TRUE)
}
