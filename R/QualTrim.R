#' Trim adapters and/or read ends across all files in a LongReadQC object
#'
#' Runs a trimming workflow for each file in a \code{LongReadQC} object. If
#' \code{AdapterSeq} is provided, adapter trimming is performed first. If
#' \code{Start} and/or \code{End} are provided, reads are then trimmed from the
#' 5' and/or 3' ends. Trimming results and output paths are recorded per file in
#' \code{qc_obj@metadata}.
#'
#' @param qc_obj A \code{LongReadQC} object.
#' @param OutDir Character(1). Output directory where trimmed read files will be
#'   written. Created if it does not exist.
#' @param AdapterSeq Character(1) or \code{NULL}. Adapter sequence to remove. If
#'   \code{NULL}, adapter trimming is skipped.
#' @param MaxMismatchEnd Integer(1). Maximum mismatches allowed when matching the adapter.
#' @param MinOverlapEnd Integer(1). Maximum distance from either end within which an
#'   adapter hit is treated as an end adapter and trimmed.
#' @param MinInternalDistance Integer(1). Minimum distance from either end for an
#'   adapter hit to be considered internal (chimera splitting).
#' @param MinFragmentLength Integer(1). Minimum fragment length to keep after trimming.
#' @param Start Integer(1) or \code{NULL}. Number of bases to remove from the 5' end.
#'   If \code{NULL}, no 5' trimming is applied.
#' @param End Integer(1) or \code{NULL}. Number of bases to remove from the 3' end.
#'   If \code{NULL}, no 3' trimming is applied.
#' @param OutFileType Character vector specifying output format(s). Supported values
#'   are \code{"fastq"}, \code{"fasta"}, and \code{"bam"}. Multiple formats may be
#'   requested.
#' @param verbose Logical(1). If \code{TRUE}, prints progress messages.
#'
#' @return A \code{LongReadQC} object with per-file trimming results stored in
#'   \code{qc_obj@metadata[[file]]} (e.g., \code{$adapter_summary}, \code{$trim_summary}).
#'
#' @seealso \code{\link{RemoveAdapter}} and \code{\link{TrimLong}}.
#'
#' @examples
#' \dontrun{
#' qc_trim <- QualTrim(
#'   qc_obj,
#'   OutDir = "trimmed_reads",
#'   AdapterSeq = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
#'   Start = 10,
#'   End = 5,
#'   OutFileType = "fastq",
#'   verbose = TRUE
#' )
#' }
#'
#' @export
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

