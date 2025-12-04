TrimLong <- function(qc_obj,
                     Start = NULL,
                     End = NULL,
                     FilePath = NULL,
                     OutDir = ".",
                     OutFile = NULL,
                     OutFileType = c("fastq")) {
  
  # Choose which file to trim
  fpath <- if (!is.null(FilePath)) FilePath else qc_obj@files[[1]]
  OriginalFPath <- qc_obj@files[[1]]
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  
  # Import reads
  reads <- ImportFile(fpath)
  
  # Compute trimming
  n <- Biostrings::width(reads)
  trim_start <- rep(if (!is.null(Start)) Start + 1 else 1, length(n))
  trim_end   <- if (!is.null(End)) n - End else n
  
  trimmed_reads <- IRanges::narrow(reads, start = trim_start, end = trim_end)
  
  # Output naming
  base_name <- if (!is.null(OutFile)) RemoveExt(OutFile) else RemoveExt(fpath)
  
  out_paths <- WriteReadOutputs(
    reads = trimmed_reads,
    base_name = paste0(base_name, "_trimmed"),
    OutDir = OutDir,
    OutFileType = OutFileType
  )
  
  # metadata storage (per file path)
  qc_obj@metadata[[OriginalFPath]]$trim_summary <- list(
    file_out = out_paths[[1]],
    reads_before = length(reads),
    reads_after = length(trimmed_reads),
    start_trim = ifelse(is.null(Start), 0, Start),
    end_trim = ifelse(is.null(End), 0, End),
    yield_after = sum(Biostrings::width(trimmed_reads)),
    trim_date = Sys.time()
  )
  
  base::return(qc_obj)
}

