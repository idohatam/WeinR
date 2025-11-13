RemoveExt <- function(x) {
  sub("\\.(fastq|fq|fasta|fa|bam)(\\.gz)?$", "", basename(x), ignore.case = TRUE)
}

WriteReadOutputs <- function(reads, base_name, OutDir = ".", OutFileType = c("fastq")) {
  # Ensure output directory exists
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  
  out_paths <- character(0)
  
  for (otype in OutFileType) {
    otype <- tolower(otype)
    if (!otype %in% c("fasta", "fastq", "bam"))
      stop("Unsupported output type: ", otype)
    
    out_ext <- switch(otype,
                      fasta = "fasta",
                      fastq = "fastq",
                      bam = "bam")
    
    out_name <- paste0(base_name, ".", out_ext)
    out_path <- file.path(OutDir, out_name)
    
    if (otype == "fasta") {
      Biostrings::writeXStringSet(as(reads, "DNAStringSet"),
                                  filepath = out_path, format = "fasta")
      
    } else if (otype == "fastq") {
      Biostrings::writeXStringSet(reads,
                                  filepath = out_path, format = "fastq",
                                  compress = grepl("\\.gz$", out_path))
      
    } else if (otype == "bam") {
      qname <- names(reads)
      seq <- as.character(reads)
      qual <- as.character(Biostrings::PhredQuality(Biostrings::quality(reads)))
      
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
      write.table(bam_df, file = temp_sam, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = FALSE, append = TRUE)
      Rsamtools::asBam(temp_sam,
                       destination = tools::file_path_sans_ext(out_path),
                       overwrite = TRUE)
      unlink(temp_sam)
    }
    
    message(length(reads), " reads written to ", out_path)
    out_paths <- c(out_paths, out_path)
  }
  
  return(out_paths)
}