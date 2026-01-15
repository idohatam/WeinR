#' Validate an input file path and infer file type (internal)
#'
#' Internal helper that checks a single file path, warns on extremely large files,
#' and returns the inferred input type based on the filename extension.
#'
#' Supported extensions are BAM (\code{.bam}), FASTQ (\code{.fastq}, \code{.fq}),
#' and gzipped FASTQ (\code{.fastq.gz}, \code{.fq.gz}).
#'
#' @param filePath Character(1). Path to an input file.
#'
#' @return A length-1 character string: one of \code{"bam"}, \code{"fastq"},
#'   or \code{"fastq.gz"}.
#'
#' @details
#' This function is intended for internal use and is not part of the public API.
#' It errors on invalid paths or unsupported extensions.
#'
#' @examples
#' # FASTQ example (creates a temporary file)
#' fq <- tempfile(fileext = ".fastq")
#' file.create(fq)
#' CheckFile(fq)
#'
#' # gzipped FASTQ example (still just a temp file name here)
#' fqgz <- tempfile(fileext = ".fastq.gz")
#' file.create(fqgz)
#' CheckFile(fqgz)
#'
#' # BAM example
#' bam <- tempfile(fileext = ".bam")
#' file.create(bam)
#' CheckFile(bam)
#'
#' @keywords internal
CheckFile <- function(filePath) 
{
  
  # check if file path is valid, if not stop with an error message
  if (is.null(filePath)) stop("File path is NULL.")
  if (!is.character(filePath)) stop("File path must be a character string.")
  if (length(filePath) != 1L) stop("File path must be a single path (length 1).")
  if (is.na(filePath)) stop("File path is NA.")
  if (!nzchar(trimws(filePath))) stop("File path is empty.")
  
  # check if file exists, if not stop with an error message
  if (!file.exists(filePath)) stop(paste("File does not exist: ", filePath))
  if (dir.exists(filePath)) stop(paste("Expected a file path, got a directory: ", filePath))
  
  # large file warning (edge case)
  sz <- file.info(filePath)$size
  if (!is.na(sz) && sz >= 20 * 1024^3) {
    warning("Large file (>=20GB). Processing may be slow.")
  }
  
  if (grepl("\\.bam$", filePath, ignore.case = TRUE)) 
  {
    return("bam")
  } 
  else if (grepl("\\.(fastq|fq)\\.gz$", filePath, ignore.case = TRUE)) 
  {
    return("fastq.gz")
  } 
  else if (grepl("\\.(fastq|fq)$", filePath, ignore.case = TRUE)) 
  {
    return("fastq")
  } 
  else 
  {
    # if file isn't supported, stop with an error message
    stop("Unsupported file type. Expected .bam or .fastq/.fq (optionally .gz).")
  }
}
