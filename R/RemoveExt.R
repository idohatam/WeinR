#' Remove sequencing file extensions from a path (internal)
#'
#' Internal helper that removes common sequencing-related file extensions
#' and sanitizes the resulting filename for safe use in output directories.
#'
#' @param x Character vector of file paths or filenames.
#'
#' @return A character vector with sequencing file extensions removed and
#' sanitized (non-alphanumeric characters replaced with "_").
#'
#' @keywords internal
#' @noRd
RemoveExt <- function(x) {
  base <- sub(
    "\\.(fastq|fq|fasta|fa|bam)(\\.gz)?$",
    "",
    basename(x),
    ignore.case = TRUE
  )
  
  # sanitize to portable filename
  base::gsub("[^A-Za-z0-9._-]+", "_", base)
}

