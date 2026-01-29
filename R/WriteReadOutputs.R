#' Write reads to one or more output formats (internal)
#'
#' Internal helper that writes a read set to disk in one or more formats.
#' Supported output formats are FASTQ, FASTA, and BAM.
#'
#' @param reads A \code{Biostrings::QualityScaledDNAStringSet} (or compatible
#'   \code{XStringSet}) containing reads to write.
#' @param base_name Character(1). Output filename stem (without extension).
#' @param OutDir Character(1). Output directory. Created if it does not exist.
#' @param OutFileType Character vector specifying output format(s). Supported values
#'   are \code{"fastq"}, \code{"fasta"}, and \code{"bam"}. Multiple formats may be
#'   requested.
#' @param Verbose Logical(1). If \code{TRUE}, print a message for each output file written.
#'
#' @return A character vector of output file paths (one per requested format).
#'
#' @details
#' FASTA/FASTQ outputs are written using \code{Biostrings::writeXStringSet()}.
#' BAM output is created by writing a temporary SAM file and converting it via
#' \code{Rsamtools::asBam()}.
#'
#' @examples
#' \dontrun{
#' reads <- Biostrings::QualityScaledDNAStringSet(
#'   Biostrings::DNAStringSet(c("ACGT", "AAAA")),
#'   Biostrings::BStringSet(c("IIII", "####"))
#' )
#' names(reads) <- c("r1", "r2")
#'
#' out <- WriteReadOutputs(
#'   reads,
#'   base_name   = "toy",
#'   OutDir      = tempdir(),
#'   OutFileType = c("fasta", "fastq"),
#'   Verbose     = TRUE
#' )
#' out
#' }
#'
#' @keywords internal
WriteReadOutputs <- function(reads,
                             base_name,
                             OutDir = ".",
                             OutFileType = c("fastq"),
                             Verbose = TRUE) {
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  
  out_paths <- character(0)
  
  for (otype in OutFileType) {
    otype <- tolower(otype)
    if (!otype %in% c("fasta", "fastq", "bam"))
      stop("Unsupported output type: ", otype)
    
    out_ext <- switch(otype,
                      fasta = "fasta",
                      fastq = "fastq",
                      bam   = "bam"
    )
    
    out_name <- paste0(base_name, ".", out_ext)
    out_path <- file.path(OutDir, out_name)
    
    if (otype == "fasta") {
      Biostrings::writeXStringSet(
        as(reads, "DNAStringSet"),
        filepath = out_path,
        format = "fasta"
      )
      
    } else if (otype == "fastq") {
      Biostrings::writeXStringSet(
        reads,
        filepath = out_path,
        format = "fastq"
      )
      
    } else if (otype == "bam") {
      qname <- names(reads)
      seq   <- as.character(reads)
      qual  <- as.character(Biostrings::PhredQuality(Biostrings::quality(reads)))
      
      bam_df <- data.frame(
        qname = qname,
        flag = 4,
        rname = "*",
        pos = 0,
        mapq = 255,
        cigar = "*",
        rnext = "*",
        pnext = 0,
        tlen = 0,
        seq = seq,
        qual = qual,
        stringsAsFactors = FALSE
      )
      
      temp_sam <- tempfile(fileext = ".sam")
      writeLines("@HD\tVN:1.6\tSO:unsorted", con = temp_sam)
      write.table(
        bam_df,
        file = temp_sam,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        append = TRUE
      )
      
      Rsamtools::asBam(
        temp_sam,
        destination = tools::file_path_sans_ext(out_path),
        overwrite = TRUE
      )
      unlink(temp_sam)
    }
    
    if (Verbose) message(length(reads), " reads written to ", out_path)
    out_paths <- c(out_paths, out_path)
  }
  
  return(out_paths)
}

