#' Remove sequencing file extensions from a path
#'
#' Internal helper that removes common sequencing-related file extensions from
#' file paths or filenames.
#'
#' Supported extensions include FASTQ, FASTA, and BAM, optionally gzipped.
#'
#' @param x Character vector of file paths or filenames.
#'
#' @return A character vector with sequencing file extensions removed.
#'
#' @details
#' The following extensions are stripped (case-insensitive):
#' \code{.fastq}, \code{.fq}, \code{.fasta}, \code{.fa}, \code{.bam},
#' and their gzipped forms (e.g., \code{.fastq.gz}).
#'
#' This function is intended for internal use and is not part of the public API.
#'
#' @examples
#' \dontrun{
#' WeinR:::RemoveExt("reads.fastq")
#' WeinR:::RemoveExt("reads.fastq.gz")
#' WeinR:::RemoveExt("/path/to/reads.FA")
#' WeinR:::RemoveExt(c("a.fq", "b.bam.gz"))
#' }
#'
#' @keywords internal
RemoveExt <- function(x) {
  sub("\\.(fastq|fq|fasta|fa|bam)(\\.gz)?$", "", basename(x), ignore.case = TRUE)
}
