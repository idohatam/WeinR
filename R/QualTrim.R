QualTrim <- function(qc_obj,
                     OutDir = ".",
                     AdapterSeq = NULL,
                     MaxMismatchEnd = 3,
                     MinOverlapEnd = 20,
                     MinInternalDistance = 100,
                     MinFragmentLength = 200,
                     Start = NULL,
                     End = NULL,
                     OutFileType = c("fastq"),
                     verbose = TRUE) {
  
  message("Starting trimming workflow for ", length(qc_obj@files), " file(s)...")
  
  # Ensure output directory exists
  OutDir <- normalizePath(OutDir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  
  # Loop through each file
  for (fpath in qc_obj@files) {
    
    message("\n→ Processing file: ", fpath)
    
    # Build single-file subobject
    single_qc <- qc_obj
    single_qc@files <- fpath
    single_qc@metrics[[fpath]] <- qc_obj@metrics[[fpath]]
    single_qc@summary_metrics <- qc_obj@summary_metrics[
      qc_obj@summary_metrics$file == fpath,
      , drop = FALSE
    ]
    
    ## Adapter trimming (RemoveAdapter)
    if (!is.null(AdapterSeq)) {
      if (verbose) message("   - Running adapter trimming...")
      
      single_qc <- RemoveAdapter(
        qc_obj = single_qc,
        adapterSeq = AdapterSeq,
        MaxMismatchEnd = MaxMismatchEnd,
        MinOverlapEnd = MinOverlapEnd,
        MinInternalDistance = MinInternalDistance,
        MinFragmentLength = MinFragmentLength,
        OutDir = OutDir,
        OutFileType = OutFileType,
        verbose = verbose
      )
    }
    
    ## Step 2: 5'/3' trimming
    if (!is.null(Start) || !is.null(End)) {
      
      if (verbose) message("   - Running start/end trimming...")
      
      # Try to extract adapter output
      adapter_out <- NULL
      
      if (!is.null(AdapterSeq)) {
        ad_sum <- single_qc@metadata[[fpath]]$adapter_summary
        
        if (!is.null(ad_sum) && !is.null(ad_sum$output_paths)) {
          # Prefer FASTQ if present
          fastq_idx <- grep("\\.fastq$", ad_sum$output_paths, ignore.case = TRUE)
          if (length(fastq_idx) > 0) {
            adapter_out <- ad_sum$output_paths[fastq_idx[1]]
          } else {
            adapter_out <- ad_sum$output_paths[1]
          }
        }
      }
      
      # Use adapter-trimmed file OR original file
      trim_input <- if (!is.null(adapter_out)) adapter_out else NULL
      
      single_qc <- TrimLong(
        qc_obj      = single_qc,
        Start       = Start,
        End         = End,
        FilePath    = trim_input,
        OutDir      = OutDir,
        OutFileType = OutFileType,
        OutFile     = NULL
      )
    }
    
    ## Store metadata
    qc_obj@metadata[[fpath]] <- single_qc@metadata[[fpath]]
  }
  
  message("\nTrimming workflow complete for all files.")
  
  return(qc_obj)
}

