#' Write reads to one or more output formats (internal)
#'
#' Internal helper that writes a read set to disk in one or more formats.
#' Supported output formats are FASTQ, FASTA, and BAM.
#'
#' @param reads An \code{XStringSet}-compatible object containing reads to write.
#'   For \code{"fastq"} and \code{"bam"} outputs, \code{reads} must be a
#'   \code{Biostrings::QualityScaledDNAStringSet}.
#' @param BaseName Character(1). Output filename stem (without extension).
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
#' FASTA is written with \code{Biostrings::writeXStringSet()}. FASTQ is written with
#' \code{ShortRead::writeFastq()} for efficient streaming. BAM output is created by
#' writing a temporary SAM file and converting it via \code{Rsamtools::asBam()}.
#'
#' @keywords internal
WriteReadOutputs <- function(reads,
                             BaseName,
                             OutDir = ".",
                             OutFileType = c("fastq"),
                             Verbose = TRUE) {
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)
  
  if (!inherits(reads, "XStringSet"))
    stop("WriteReadOutputs(): reads must be an XStringSet-compatible object.", call. = FALSE)
  
  n_reads <- length(reads)
  if (n_reads == 0L) return(character(0))
  
  # Ensure IDs exist and are safe for headers/QNAME
  ids <- names(reads)
  if (is.null(ids) || length(ids) != n_reads || anyNA(ids) || any(!nzchar(ids))) {
    ids <- paste0("read_", seq_len(n_reads))
  }
  ids <- gsub("[\\r\\n\\t]+", "_", ids)
  
  out_paths <- character(0)
  
  for (otype in OutFileType) {
    otype <- tolower(otype)
    if (!otype %in% c("fasta", "fastq", "bam"))
      stop("WriteReadOutputs(): unsupported output type: ", otype, call. = FALSE)
    
    out_ext  <- switch(otype, fasta = "fasta", fastq = "fastq", bam = "bam")
    out_path <- file.path(OutDir, paste0(BaseName, ".", out_ext))
    
    if (otype == "fasta") {
      
      Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(reads),
        filepath = out_path,
        format   = "fasta"
      )
      
    } else if (otype == "fastq") {
      
      if (!inherits(reads, "QualityScaledDNAStringSet")) {
        stop("WriteReadOutputs(): FASTQ output requires a QualityScaledDNAStringSet.", call. = FALSE)
      }
      
      sr <- ShortRead::ShortReadQ(
        sread   = Biostrings::DNAStringSet(reads),
        quality = ShortRead::FastqQuality(Biostrings::quality(reads)),
        id      = Biostrings::BStringSet(ids)
      )
      
      ShortRead::writeFastq(sr, out_path, compress = FALSE)
      
    } else if (otype == "bam") {
      
      if (!inherits(reads, "QualityScaledDNAStringSet")) {
        stop("WriteReadOutputs(): BAM output requires a QualityScaledDNAStringSet.", call. = FALSE)
      }
      
      # QNAME must not contain whitespace: keep only first token
      qname <- sub("\\s.*$", "", ids)
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
      
      writeLines(c(
        "@HD\tVN:1.6\tSO:unsorted",
        "@SQ\tSN:*\tLN:0"
      ), con = temp_sam)
      
      utils::write.table(
        bam_df, file = temp_sam, sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = FALSE, append = TRUE
      )
      
      Rsamtools::asBam(
        temp_sam,
        destination = tools::file_path_sans_ext(out_path),
        overwrite   = TRUE
      )
      unlink(temp_sam)
    }
    
    if (isTRUE(Verbose)) message(n_reads, " reads written to ", out_path)
    out_paths <- c(out_paths, out_path)
  }
  
  out_paths
}
