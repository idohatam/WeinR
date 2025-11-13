TrimLong <- function(qc_obj,
                     Start = NULL,
                     End = NULL,
                     OutDir = ".",
                     OutFileType = c("fastq"),
                     OutFile = NULL) {
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  
  # Correct file info extraction
  fpath <- qc_obj@files[[1]]
  fname <- basename(fpath)   # Now only for printing
  
  message("Trimming reads for file: ", fname)
  
  # Import reads
  reads <- ImportFile(fpath)
  
  # Compute trimming boundaries
  n <- Biostrings::width(reads)
  trim_start <- rep(if (!is.null(Start)) Start + 1 else 1, length(n))
  trim_end   <- if (!is.null(End)) n - End else n
  
  # Apply trimming
  trimmed_reads <- IRanges::narrow(reads, start = trim_start, end = trim_end)
  
  # Output naming (fixed)
  base_name <- if (!is.null(OutFile)) {
    RemoveExt(OutFile)
  } else {
    RemoveExt(fpath)   # use full path consistently
  }
  
  out_paths <- WriteReadOutputs(
    reads = trimmed_reads,
    base_name = paste0(base_name, "_trimmed"),
    OutDir = OutDir,
    OutFileType = OutFileType
  )
  
  # Metadata stored using the file path as the key
  qc_obj@metadata[[fpath]]$trim_summary <- list(
    start_trimmed = ifelse(is.null(Start), 0, Start),
    end_trimmed   = ifelse(is.null(End), 0, End),
    total_reads   = length(trimmed_reads),
    avg_length    = mean(Biostrings::width(trimmed_reads)),
    yield_after   = sum(Biostrings::width(trimmed_reads)),
    output_paths  = out_paths,
    trim_date     = Sys.time()
  )
  
  message(
    "Trimming complete for ", fname, ": ",
    qc_obj@metadata[[fpath]]$trim_summary$total_reads, " reads retained. ",
    "(start=", ifelse(is.null(Start), 0, Start),
    ", end=", ifelse(is.null(End), 0, End), ")"
  )
  
  base::return(qc_obj)
}
