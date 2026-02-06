#' Write reads to one or more output formats (internal)
#'
#' Internal helper that writes a read set to disk in one or more formats.
#' Supported output formats are FASTQ, FASTA, and BAM.
#'
#' @param reads An \code{XStringSet}-compatible object containing reads to write.
#'   For \code{"fastq"} and \code{"bam"} outputs, \code{reads} must be a
#'   \code{Biostrings::QualityScaledDNAStringSet}.
#' @param base_name Character(1). Output filename stem (without extension).
#' @param OutDir Character(1). Output directory. Created if it does not exist.
#' @param OutFileType Character vector specifying output format(s). Supported values
#'   are \code{"fastq"}, \code{"fasta"}, and \code{"bam"}. Multiple formats may be
#'   requested.
#' @param Verbose Logical(1). If \code{TRUE}, print a message for each output file written.
#'
#' @return A character vector of output file paths (one per requested format).
#'   Returns \code{character(0)} if \code{reads} contains \code{0} sequences.
#'
#' @details
#' FASTA and FASTQ outputs are written using \code{writeLines()} with simple
#' record formatting. BAM output is created by writing a temporary SAM file and
#' converting it via \code{Rsamtools::asBam()}.
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
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)
  
  if (!inherits(reads, "XStringSet"))
    stop("WriteReadOutputs(): reads must be an XStringSet-compatible object.")
  
  n_reads <- length(reads)
  if (n_reads == 0L) return(character(0))
  
  out_paths <- character(0)
  
  # Ensure names exist (FASTA/FASTQ/BAM require IDs)
  ids <- names(reads)
  if (is.null(ids) || any(is.na(ids)) || any(!nzchar(ids))) {
    ids <- paste0("read_", seq_len(n_reads))
  }
  
  for (otype in OutFileType) {
    otype <- tolower(otype)
    if (!otype %in% c("fasta", "fastq", "bam"))
      stop("WriteReadOutputs(): unsupported output type: ", otype)
    
    out_ext <- switch(otype,
                      fasta = "fasta",
                      fastq = "fastq",
                      bam   = "bam")
    
    out_path <- file.path(OutDir, paste0(base_name, ".", out_ext))
    
    # FASTA (writeLines)
    if (otype == "fasta") {
      seqs <- as.character(Biostrings::DNAStringSet(reads))
      fasta_lines <- as.vector(rbind(paste0(">", ids), seqs))
      writeLines(fasta_lines, con = out_path)
      
    } else if (otype == "fastq") {
      # FASTQ requires qualities
      if (!inherits(reads, "QualityScaledDNAStringSet")) {
        stop("WriteReadOutputs(): FASTQ output requires a QualityScaledDNAStringSet.")
      }
      
      seqs  <- as.character(Biostrings::DNAStringSet(reads))
      quals <- as.character(Biostrings::quality(reads))
      
      fastq_lines <- character(4L * n_reads)
      idx <- seq.int(1L, by = 4L, length.out = n_reads)
      fastq_lines[idx]     <- paste0("@", ids)
      fastq_lines[idx + 1] <- seqs
      fastq_lines[idx + 2] <- "+"
      fastq_lines[idx + 3] <- quals
      
      writeLines(fastq_lines, con = out_path)
      
    } else if (otype == "bam") {
      # BAM requires qualities
      if (!inherits(reads, "QualityScaledDNAStringSet")) {
        stop("WriteReadOutputs(): BAM output requires a QualityScaledDNAStringSet.", call. = FALSE)
      }
      
      n <- length(reads)
      qname <- names(reads)
      if (is.null(qname) || any(is.na(qname)) || any(!nzchar(qname))) {
        qname <- paste0("read_", seq_len(n))
      }
      
      # QNAME must not contain whitespace: keep only first token
      qname <- sub("\\s.*$", "", qname)
      qname <- gsub("[\\t\\r\\n ]+", "_", qname)
      
      seq  <- as.character(Biostrings::DNAStringSet(reads))
      qual <- as.character(Biostrings::quality(reads))
      
      bam_df <- data.frame(
        qname = qname,
        flag  = 4L,
        rname = "*",
        pos   = 0L,
        mapq  = 0L,
        cigar = "*",
        rnext = "*",
        pnext = 0L,
        tlen  = 0L,
        seq   = seq,
        qual  = qual,
        stringsAsFactors = FALSE
      )
      
      temp_sam <- tempfile(fileext = ".sam")
      
      # Minimal SAM header
      writeLines(c(
        "@HD\tVN:1.6\tSO:unsorted",
        "@SQ\tSN:*\tLN:0"
      ), con = temp_sam)
      
      write.table(bam_df, file = temp_sam, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = FALSE, append = TRUE)
      
      Rsamtools::asBam(
        temp_sam,
        destination = tools::file_path_sans_ext(out_path),
        overwrite = TRUE
      )
      unlink(temp_sam)
    }
    
    if (isTRUE(Verbose)) message(n_reads, " reads written to ", out_path)
    out_paths <- c(out_paths, out_path)
  }
  
  out_paths
}

