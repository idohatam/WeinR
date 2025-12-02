#' Run filtering workflow across all files in a LongReadQC object
#'
#' @param qc_obj A LongReadQC object (with metrics already filled by QualMat()).
#' @param OutDir Directory where filtered read files will be saved.
#' @param MinAvgQS Minimum per-read average quality score threshold.
#' @param MinLength Minimum read length threshold.
#' @param MaxNumberNs Maximum allowed number of Ns per read.
#' @param OutFileType Output format(s): "fastq", "fasta", or "bam".
#'
#' @return The updated LongReadQC object with filtering results in metadata.
#' @export

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
    
    # Build single-file subobject
    single_qc <- qc_obj
    single_qc@files <- fpath
    single_qc@metrics[[fpath]] <- qc_obj@metrics[[fpath]]
    single_qc@summary_metrics <- qc_obj@summary_metrics[
      qc_obj@summary_metrics$file == fpath,
      , drop = FALSE
    ]
    
    # Run filtering
    filtered_qc <- FilterLong(
      qc_obj = single_qc,
      MinAvgQS = MinAvgQS,
      MinLength = MinLength,
      MaxNumberNs = MaxNumberNs,
      OutFileType = OutFileType,
      OutDir = OutDir
    )
    
    # If filtering produced no metadata, skip
    if (is.null(filtered_qc@metadata[[fpath]]$filter_summary)) {
      message("No filtering summary for: ", fpath)
      next
    }
    
    # Save ONLY this file’s metadata
    qc_obj@metadata[[fpath]] <- filtered_qc@metadata[[fpath]]
  }
  
  message("\nFiltering complete for all files.")
  
  return(qc_obj)
}
