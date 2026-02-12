#' Validate an input file path and infer file type
#'
#' Internal helper that validates a single input file path, issues a warning
#' for large files, and infers the input type from the filename
#' extension.
#'
#' Supported extensions are BAM (\code{.bam}), FASTQ (\code{.fastq}, \code{.fq}),
#' and gzipped FASTQ (\code{.fastq.gz}, \code{.fq.gz}).
#'
#' @param filePath Character scalar. Path to an input file.
#'
#' @return A length-1 character string: one of \code{"bam"},
#'   \code{"fastq"}, or \code{"fastq.gz"}.
#'
#' @details
#' This function is intended for internal use only and is not part of the
#' public API. It errors on invalid paths or unsupported file extensions.
#'
#' @keywords internal
CheckFile <- function(filePath) 
{
  
  # check if file path is valid, if not stop with an error message
  if (is.null(filePath)) stop("File path is NULL.", call. = FALSE)
  if (!is.character(filePath)) stop("File path must be a character string.", call. = FALSE)
  if (length(filePath) != 1L) stop("File path must be a single path (length 1).", call. = FALSE)
  if (is.na(filePath)) stop("File path is NA.", call. = FALSE)
  if (!nzchar(trimws(filePath))) stop("File path is empty.", call. = FALSE)
  
  # check if file exists, if not stop with an error message
  if (!file.exists(filePath)) stop(paste("File does not exist: ", filePath), call. = FALSE)
  if (dir.exists(filePath)) stop(paste("Expected a file path, got a directory: ", filePath), call. = FALSE)
  
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
    stop("Unsupported file type. Expected .bam or .fastq/.fq (optionally .gz).", call. = FALSE)
  }
}
