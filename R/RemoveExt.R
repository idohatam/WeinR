#' Remove sequencing file extensions from a path (internal)
#'
#' Internal helper that removes common sequencing file extensions from a file
#' path or filename. Supported extensions include FASTQ, FASTA, and BAM,
#' optionally gzipped.
#'
#' @param x Character vector of file paths or filenames.
#'
#' @return A character vector with sequencing file extensions removed.
#'
#' @details
#' The function strips the following extensions (case-insensitive):
#' \code{.fastq}, \code{.fq}, \code{.fasta}, \code{.fa}, \code{.bam},
#' and their gzipped forms (e.g., \code{.fastq.gz}).
#'
#' @examples
#' RemoveExt("reads.fastq")
#' RemoveExt("reads.fastq.gz")
#' RemoveExt("/path/to/reads.FA")
#' RemoveExt(c("a.fq", "b.bam.gz"))
#'
#' @keywords internal
RemoveExt <- function(x) {
  sub("\\.(fastq|fq|fasta|fa|bam)(\\.gz)?$", "", basename(x), ignore.case = TRUE)
}
