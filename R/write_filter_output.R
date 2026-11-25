

write_filter_output <- function(processed_chunk, base_dir, filename, timestamp, output_filetype) {
  
  n_reads <- length(processed_chunk)
  ids <- names(processed_chunk)
  if (is.null(ids) || length(ids) != n_reads || anyNA(ids) || any(!nzchar(ids))) {
    ids <- paste0("read_", seq_len(n_reads))
  }
  ids <- gsub("[\\r\\n\\t]+", "_", ids)
  
  filename_no_extension <- tools::file_path_sans_ext(filename)
  output_file <- paste0(base_dir, "/Processed_Files/", filename_no_extension, "_", timestamp)
  
  if(output_filetype == "fastq") {
    sr <- ShortRead::ShortReadQ(
      sread = Biostrings::DNAStringSet(processed_chunk),
      quality = ShortRead::FastqQuality(Biostrings::quality(processed_chunk)),
      id = Biostrings::BStringSet(ids)
    )
    ShortRead::writeFastq(sr, paste0(output_file, ".fastq"), mode = "a", compress = FALSE)
  } else if (output_filetype == "fasta") {
    Biostrings::writeXStringSet(
      Biostrings::DNAStringSet(processed_chunk),
      filepath = paste0(output_file, ".fasta"),
      format = "fasta",
      append = TRUE
    )
  } else if(output_filetype == "bam") {
    # QNAME must not contain whitespace: keep only first token
    qname <- sub("\\s.*$", "", ids)
    qname <- gsub("[\\t\\r\\n ]+", "_", qname)
    seq  <- as.character(Biostrings::DNAStringSet(processed_chunk))
    qual <- as.character(Biostrings::quality(processed_chunk))
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
    temp_sam <- paste0(output_file, ".sam")
    utils::write.table(
      bam_df, file = temp_sam, sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = FALSE, append = TRUE
    )
  } else {
    stop("WriteFilterOutput(): unsupported output type: ", output_filetype, call. = FALSE)
  }
  output_file
}