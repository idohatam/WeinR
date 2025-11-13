FilterLong <- function(qc_obj,
                       MinAvgQS = 20,
                       MinLength = 100,
                       MaxNumberNs = 2,
                       OutFileType = c("fastq"),
                       OutDir = ".") {
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  
  ## Extract file path and label
  fpath <- qc_obj@files[[1]]
  fname <- basename(fpath)
  
  ## Retrieve metrics for this file using full path
  file_metrics <- qc_obj@metrics[[fpath]]
  if (is.null(file_metrics))
    stop("Metrics not found for file: ", fpath,
         "\nMake sure QualMat() used the full path as the key.")
  
  ## Extract per-read metrics
  lengths  <- file_metrics$readLengths
  q_scores <- unlist(file_metrics$meanprQscore)
  n_counts <- unlist(file_metrics$Ncount)
  

  n <- length(lengths)
  if (length(q_scores) != n || length(n_counts) != n)
    stop("Mismatch in metric lengths — QualMat() produced inconsistent per-read lists.")
  

  keep <- (q_scores >= MinAvgQS) &
    (lengths  >= MinLength) &
    (n_counts <= MaxNumberNs)
  
  ## If nothing passes, return qc_obj unchanged
  if (!any(keep)) {
    message("No reads passed filters for: ", fname)
    base::return(qc_obj)
  }
  
  ## Load reads
  reads <- ImportFile(fpath)
  
  ## Filter reads BY INDEX (Check with Ido + Myles if names of reads should be saved in the quality matrix)
  filtered_reads <- reads[keep]
  
  ## Write output
  base_name <- paste0(RemoveExt(fpath), "_filtered")
  
  out_paths <- WriteReadOutputs(
    reads = filtered_reads,
    base_name = base_name,
    OutDir = OutDir,
    OutFileType = OutFileType
  )
  
  ## Update metadata
  qc_obj@metadata$filter_summary <- list(
    reads_before = length(reads),
    reads_after  = length(filtered_reads),
    yield_before = sum(Biostrings::width(reads)),
    yield_after  = sum(Biostrings::width(filtered_reads)),
    filter_params = list(
      MinAvgQS = MinAvgQS,
      MinLength = MinLength,
      MaxNumberNs = MaxNumberNs
    ),
    output_paths = out_paths,
    filter_date = Sys.time()
  )
  
  ## Return updated object
  base::return(qc_obj)
}

