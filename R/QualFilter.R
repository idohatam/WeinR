QualFilter <- function(qc_obj,
                       OutDir = ".",
                       MinAvgQS = 20,
                       MinLength = 100,
                       MaxNumberNs = 2,
                       OutFileType = c("fastq")) {
  
  message("Starting filtering workflow for ", length(qc_obj@files), " file(s)...")
  
  OutDir <- normalizePath(OutDir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  
  for (fpath in qc_obj@files) {
    
    message("\n→ Processing file: ", fpath)
    
    # Construct a single-file subobject
    single_qc <- qc_obj
    single_qc@files <- fpath
    
    # extract metrics by path
    single_qc@metrics[[fpath]] <- qc_obj@metrics[[fpath]]
    # subset summary_metrics by *path*
    single_qc@summary_metrics <- qc_obj@summary_metrics[
      qc_obj@summary_metrics$file == fpath,
      , drop = FALSE
    ]
    
    # run filtering
    filtered_qc <- FilterLong(
      qc_obj = single_qc,
      MinAvgQS = MinAvgQS,
      MinLength = MinLength,
      MaxNumberNs = MaxNumberNs,
      OutFileType = OutFileType,
      OutDir = OutDir
    )
    
    # no filter_summary produced (no reads passed)
    if (is.null(filtered_qc@metadata$filter_summary)) {
      message("No filtering summary for: ", fpath, " — skipping summary update.")
      next
    }
    
    # update summary metrics
   # qc_obj@summary_metrics[
      #qc_obj@summary_metrics$file == fpath,
      "filtered_yield"
    #] <- filtered_qc@metadata$filter_summary$yield_after
    
    # Store metadata under full path
    qc_obj@metadata[[fpath]] <- filtered_qc@metadata
  }
  
  message("\nFiltering complete for all files.")
  
  return(qc_obj)
}
