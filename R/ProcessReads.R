#' Process long-read files through filtering, adapter removal, end trimming, and reporting
#'
#' Runs an end-to-end processing workflow over all files in a `LongReadQC` object.
#' For each input file, the workflow can (optionally) apply read filtering
#' (`FilterLong()`), adapter trimming (`RemoveAdapter()`), and 5'/3' end trimming
#' (`TrimLong()`), chaining steps via temporary FASTQ files when intermediate writing
#' is enabled. After processing all files, an HTML report can be generated via
#' `Reporter()`.
#'
#' @param qc_obj A `LongReadQC` object containing one or more input files in `qc_obj@files`,
#'   with corresponding metrics already computed in `qc_obj@metrics` and
#'   `qc_obj@summary_metrics`.
#' @param filter Logical(1). If `TRUE`, run `FilterLong()` prior to other steps.
#' @param MinAvgQS Numeric(1). Minimum mean per-read quality score to keep a read
#'   (passed to `FilterLong()`).
#' @param MinLength Integer(1). Minimum read length to keep a read (passed to `FilterLong()`).
#' @param MaxNumberNs Integer(1). Maximum number of `N` bases allowed per read
#'   (passed to `FilterLong()`).
#' @param AdapterSeq Character(1) or `NULL`. If not `NULL`, run `RemoveAdapter()` using
#'   this adapter sequence.
#' @param MaxMismatchEnd Integer(1). Maximum mismatches allowed in end-adapter matches
#'   (passed to `RemoveAdapter()`).
#' @param MinOverlapEnd Integer(1). Minimum overlap length for end-adapter matches
#'   (passed to `RemoveAdapter()`).
#' @param MinInternalDistance Integer(1). Minimum distance between internal adapter sites
#'   (passed to `RemoveAdapter()`).
#' @param MinFragmentLength Integer(1). Minimum fragment length to keep after internal
#'   adapter trimming (passed to `RemoveAdapter()`).
#' @param Start Integer(1) or `NULL`. If not `NULL`, trim this many bases from the 5' end
#'   (passed to `TrimLong()`).
#' @param End Integer(1) or `NULL`. If not `NULL`, trim this many bases from the 3' end
#'   (passed to `TrimLong()`).
#' @param OutFileType Character vector specifying output format(s) to write when outputs
#'   are written (e.g., \code{"fastq"}). Supported values are \code{"fastq"},
#'   \code{"bam"}, and \code{"fasta"} (as implemented by \code{WriteReadOutputs()}).
#' @param outpath Character(1). Basename for the report output (without extension).
#'   Reports are written to `file.path(getwd(), "WeinR_Outputs", "Reports", outpath)`,
#'   producing `.Rmd` and `.html`.
#' @param title Character(1). Title used in the generated report. Defaults to `outpath`.
#' @param render_report Logical(1). If \code{TRUE}, generate an HTML report with \code{Reporter()}.
#' @param force Logical(1). If `FALSE` and an existing report `.Rmd` or `.html` already exists
#'   for `outpath`, the function errors. If `TRUE`, allows overwrite.
#' @param verbose Logical(1). If `TRUE`, print progress messages for each file and step.
#' @param KeepIntermediates Logical(1). If `TRUE`, keep intermediate FASTQ files produced
#'   between steps (when step functions are called with `WriteIntermediate = TRUE`).
#'   If `FALSE`, intermediates are cleaned up after each file.
#' @param FinalSuffix Character(1). Suffix used for the final output file produced by the
#'   last step that runs (e.g., `"processed"`).
#'
#' @return The updated `LongReadQC` object with per-file metadata updated in `qc_obj@metadata`.
#'   Processed reads are written under `file.path(getwd(), "WeinR_Outputs", "Processed_Files")`
#'   by the underlying step functions.
#'
#' @details
#' **Workflow order (per file):**
#' \enumerate{
#'   \item Optional filtering with `FilterLong()` if `filter = TRUE`.
#'   \item Optional adapter trimming with `RemoveAdapter()` if `AdapterSeq` is provided.
#'   \item Optional 5'/3' trimming with `TrimLong()` if `Start` and/or `End` are provided.
#' }
#'
#' Steps may be chained using temporary FASTQ files stored as
#' `attr(single_qc, "._tmp_fastq")` produced by `FilterLong()` / `RemoveAdapter()`. 
#' Temporary files are cleaned up automatically after each file
#'
#' If a step results in `reads_after == 0L` (as recorded in the step summary in metadata),
#' remaining steps for that file are skipped and the workflow continues to the next input file.
#'
#' @seealso [FilterLong()], [RemoveAdapter()], [TrimLong()], [Reporter()]
#'
#' @examples
#' \dontrun{
#' # qc_obj is assumed to be a LongReadQC object with files + metrics already populated
#'
#' # 1) Filter only + report
#' qc2 <- ProcessReads(qc_obj,
#'   filter = TRUE,
#'   MinAvgQS = 20,
#'   MinLength = 100,
#'   MaxNumberNs = 2,
#'   outpath = "qc_filtered",
#'   render_report = TRUE
#' )
#'
#' # 2) Filter + adapter trimming + end trimming (no report)
#' qc3 <- ProcessReads(qc_obj,
#'   filter = TRUE,
#'   AdapterSeq = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
#'   Start = 10,
#'   End = 10,
#'   outpath = "qc_full",
#'   render_report = FALSE
#' )
#' }
#'
#' @export
ProcessReads <- function(qc_obj,
                         filter = FALSE,
                         MinAvgQS = 20,
                         MinLength = 200,
                         MaxNumberNs = 2,
                         AdapterSeq = NULL,
                         MaxMismatchEnd = 3,
                         MinOverlapEnd = 20,
                         MinInternalDistance = 100,
                         MinFragmentLength = 200,
                         Start = NULL,
                         End = NULL,
                         OutFileType = c("fastq"),
                         outpath,
                         title = outpath,
                         render_report = TRUE,
                         force = FALSE,
                         verbose = TRUE,
                         KeepIntermediates = FALSE,
                         FinalSuffix = "processed") {
  
  message("Starting processing workflow for ", length(qc_obj@files), " file(s)...")
  
  base_dir <- file.path(getwd(), "WeinR_Outputs")
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  
  processed_dir <- file.path(base_dir, "Processed_Files")
  processed_dir <- normalizePath(processed_dir, winslash = "/", mustWork = FALSE)
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
  
  # where to put intermediates if not keeping them
  tmp_dir <- file.path(tempdir(), "WeinR_intermediate")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (fpath in qc_obj@files) {
    
    key <- basename(fpath)
    
    tryCatch({
      
      if (verbose) message("\n-> Processing file: ", fpath)
      
      single_qc <- qc_obj
      single_qc@files <- fpath
      single_qc@metrics[[key]] <- qc_obj@metrics[[key]]
      single_qc@summary_metrics <- qc_obj@summary_metrics[
        qc_obj@summary_metrics$file == key, , drop = FALSE
      ]
      
      tmp1 <- NULL
      tmp2 <- NULL
      
      on.exit({
        if (!is.null(tmp1) && file.exists(tmp1)) unlink(tmp1)
        if (!is.null(tmp2) && file.exists(tmp2)) unlink(tmp2)
      }, add = TRUE)
      
      current_input <- fpath
      
      # Decide if there WILL be a later step (used to choose output dir + suffix)
      will_run_adapter <- !is.null(AdapterSeq)
      will_run_trim    <- !is.null(Start) || !is.null(End)
      
      
      
      # 1. Filter
      if (isTRUE(filter)) {
        if (verbose) message("   - Filtering...")
        
        is_final <- !(will_run_adapter || will_run_trim)
        single_qc <- FilterLong(
          qc_obj = single_qc,
          MinAvgQS = MinAvgQS,
          MinLength = MinLength,
          MaxNumberNs = MaxNumberNs,
          OutFileType = OutFileType,
          OutDir = processed_dir,
          WriteIntermediate = TRUE,
          KeepIntermediates = KeepIntermediates || is_final,
          OutSuffix = if (is_final) FinalSuffix else step_suffix("filtered")
        )
        
        filt_sum <- single_qc@metadata[[key]]$filter_summary
        if (!is.null(filt_sum) && isTRUE(filt_sum$reads_after == 0L)) {
          if (verbose) message("   - Filtering kept 0 reads for ", key, ". Skipping remaining steps.")
          qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
          return(NULL)
        }
        
        tmp1 <- attr(single_qc, "._tmp_fastq")
        if (!is.null(tmp1) && file.exists(tmp1)) {
          current_input <- tmp1
        } else if (!is_final) {
          if (verbose) message("   - No intermediate FASTQ from filtering for ", key, ". Skipping remaining steps.")
          qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
          return(NULL)
        }
      }
      

      # 2. Adapter Trim
      if (!is.null(AdapterSeq)) {
        if (verbose) message("   - Adapter trimming...")
        
        is_final <- !will_run_trim
        
        single_qc <- RemoveAdapter(
          qc_obj = single_qc,
          adapterSeq = AdapterSeq,
          MaxMismatchEnd = MaxMismatchEnd,
          MinOverlapEnd = MinOverlapEnd,
          MinInternalDistance = MinInternalDistance,
          MinFragmentLength = MinFragmentLength,
          FilePath = current_input,
          OutDir = processed_dir,
          OutFileType = OutFileType,
          verbose = verbose,
          WriteIntermediate = TRUE,
          KeepIntermediates = KeepIntermediates || is_final,
          OutSuffix = if (is_final) FinalSuffix else step_suffix("adaptertrimmed")
        )
        
        adap_sum <- single_qc@metadata[[key]]$adapter_summary
        if (!is.null(adap_sum) && isTRUE(adap_sum$reads_after == 0L)) {
          if (verbose) message("   - Adapter trimming kept 0 reads for ", key, ". Skipping remaining steps.")
          qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
          return(NULL)
        }
        
        tmp2 <- attr(single_qc, "._tmp_fastq")
        if (!is.null(tmp2) && file.exists(tmp2)) {
          current_input <- tmp2
        }
      }
      
      # 3. Start/End Trim (final)
      if (!is.null(Start) || !is.null(End)) {
        if (verbose) message("   - Start/end trimming...")
        
        single_qc <- TrimLong(
          qc_obj = single_qc,
          Start = Start,
          End = End,
          FilePath = current_input,
          OutDir = processed_dir,                 
          OutFileType = OutFileType,
          OutSuffix = FinalSuffix                
        )
        
        trim_sum <- single_qc@metadata[[key]]$trim_summary
        if (!is.null(trim_sum) && isTRUE(trim_sum$reads_after == 0L)) {
          if (verbose) message("   - End trimming kept 0 reads for ", key, ".")
          qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
          return(NULL)
        }
      }
      
      qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
      NULL
      
    }, error = function(e) {
      warning("Processing failed for file: ", key, "\nReason: ", conditionMessage(e))
      NULL
    })
  }
  
  message("\nProcessing complete for all files.")
  
  # 4. Report 
  reports_dir <- file.path(base_dir, "Reports")
  dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(reports_dir, outpath)
  
  rmd_path  <- paste0(report_path, ".Rmd")
  html_path <- paste0(report_path, ".html")
  
  if (!force && (file.exists(rmd_path) || file.exists(html_path))) {
    stop(
      "Report already exists at:\n",
      rmd_path, "\n",
      "Use force = TRUE to overwrite.",
      call. = FALSE
    )
  }
  
  if (isTRUE(render_report)) {
    Reporter(
      mfa         = qc_obj,
      path        = report_path,
      title       = title,
      overwrite   = TRUE,
      render_html = TRUE,
      metadata    = TRUE
    )
  }
  
  base::return(qc_obj)
}



#' Process long-read files through filtering, adapter removal, end trimming, and reporting
#'
#' Runs an end-to-end processing workflow over all files in a `LongReadQC` object.
#' For each input file, the workflow can (optionally) apply read filtering
#' (`FilterLong()`), adapter trimming (`RemoveAdapter()`), and 5'/3' end trimming
#' (`TrimLong()`), chaining steps via temporary FASTQ files when intermediate writing
#' is enabled. After processing all files, an HTML report can be generated via
#' `Reporter()`.
#'
#' @param qc_obj A `LongReadQC` object containing one or more input files in `qc_obj@files`,
#'   with corresponding metrics already computed in `qc_obj@metrics` and
#'   `qc_obj@summary_metrics`.
#' @param filter Logical(1). If `TRUE`, run `FilterLong()` prior to other steps.
#' @param MinAvgQS Numeric(1). Minimum mean per-read quality score to keep a read
#'   (passed to `FilterLong()`).
#' @param MinLength Integer(1). Minimum read length to keep a read (passed to `FilterLong()`).
#' @param MaxNumberNs Integer(1). Maximum number of `N` bases allowed per read
#'   (passed to `FilterLong()`).
#' @param AdapterSeq Character(1) or `NULL`. If not `NULL`, run `RemoveAdapter()` using
#'   this adapter sequence.
#' @param MaxMismatchEnd Integer(1). Maximum mismatches allowed in end-adapter matches
#'   (passed to `RemoveAdapter()`).
#' @param MinOverlapEnd Integer(1). Minimum overlap length for end-adapter matches
#'   (passed to `RemoveAdapter()`).
#' @param MinInternalDistance Integer(1). Minimum distance between internal adapter sites
#'   (passed to `RemoveAdapter()`).
#' @param MinFragmentLength Integer(1). Minimum fragment length to keep after internal
#'   adapter trimming (passed to `RemoveAdapter()`).
#' @param Start Integer(1) or `NULL`. If not `NULL`, trim this many bases from the 5' end
#'   (passed to `TrimLong()`).
#' @param End Integer(1) or `NULL`. If not `NULL`, trim this many bases from the 3' end
#'   (passed to `TrimLong()`).
#' @param OutFileType Character vector specifying output format(s) to write when outputs
#'   are written (e.g., \code{"fastq"}). Supported values are \code{"fastq"},
#'   \code{"bam"}, and \code{"fasta"} (as implemented by \code{WriteReadOutputs()}).
#' @param outpath Character(1). Basename for the report output (without extension).
#'   Reports are written to `file.path(getwd(), "WeinR_Outputs", "Reports", outpath)`,
#'   producing `.Rmd` and `.html`.
#' @param title Character(1). Title used in the generated report. Defaults to `outpath`.
#' @param render_report Logical(1). If \code{TRUE}, generate an HTML report with \code{Reporter()}.
#' @param force Logical(1). If `FALSE` and an existing report `.Rmd` or `.html` already exists
#'   for `outpath`, the function errors. If `TRUE`, allows overwrite.
#' @param verbose Logical(1). If `TRUE`, print progress messages for each file and step.
#' @param KeepIntermediates Logical(1). If `TRUE`, keep intermediate FASTQ files produced
#'   between steps (when step functions are called with `WriteIntermediate = TRUE`).
#'   If `FALSE`, intermediates are cleaned up after each file.
#' @param FinalSuffix Character(1). Suffix used for the final output file produced by the
#'   last step that runs (e.g., `"processed"`).
#'
#' @return The updated `LongReadQC` object with per-file metadata updated in `qc_obj@metadata`.
#'   Processed reads are written under `file.path(getwd(), "WeinR_Outputs", "Processed_Files")`
#'   by the underlying step functions.
#'
#' @details
#' **Workflow order (per file):**
#' \enumerate{
#'   \item Optional filtering with `FilterLong()` if `filter = TRUE`.
#'   \item Optional adapter trimming with `RemoveAdapter()` if `AdapterSeq` is provided.
#'   \item Optional 5'/3' trimming with `TrimLong()` if `Start` and/or `End` are provided.
#' }
#'
#' Steps may be chained using temporary FASTQ files stored as
#' `attr(single_qc, "._tmp_fastq")` produced by `FilterLong()` / `RemoveAdapter()`. 
#' Temporary files are cleaned up automatically after each file
#'
#' If a step results in `reads_after == 0L` (as recorded in the step summary in metadata),
#' remaining steps for that file are skipped and the workflow continues to the next input file.
#'
#' @seealso [FilterLong()], [RemoveAdapter()], [TrimLong()], [Reporter()]
#'
#' @examples
#' \dontrun{
#' # qc_obj is assumed to be a LongReadQC object with files + metrics already populated
#'
#' # 1) Filter only + report
#' qc2 <- ProcessReads(qc_obj,
#'   filter = TRUE,
#'   MinAvgQS = 20,
#'   MinLength = 100,
#'   MaxNumberNs = 2,
#'   outpath = "qc_filtered",
#'   render_report = TRUE
#' )
#'
#' # 2) Filter + adapter trimming + end trimming (no report)
#' qc3 <- ProcessReads(qc_obj,
#'   filter = TRUE,
#'   AdapterSeq = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
#'   Start = 10,
#'   End = 10,
#'   outpath = "qc_full",
#'   render_report = FALSE
#' )
#' }
#'
#' @export
ProcessReads <- function(qc_obj,
                         filter = FALSE,
                         MinAvgQS = 20,
                         MinLength = 200,
                         MaxNumberNs = 2,
                         AdapterSeq = NULL,
                         MaxMismatchEnd = 3,
                         MinOverlapEnd = 20,
                         MinInternalDistance = 100,
                         MinFragmentLength = 200,
                         Start = NULL,
                         End = NULL,
                         OutFileType = c("fastq"),
                         outpath,
                         title = outpath,
                         render_report = TRUE,
                         force = FALSE,
                         verbose = TRUE,
                         KeepIntermediates = FALSE,
                         FinalSuffix = "processed") {
  
  message("Starting processing workflow for ", length(qc_obj@files), " file(s)...")
  
  base_dir <- file.path(getwd(), "WeinR_Outputs")
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  
  processed_dir <- file.path(base_dir, "Processed_Files")
  processed_dir <- normalizePath(processed_dir, winslash = "/", mustWork = FALSE)
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
  
  tmp_dir <- file.path(tempdir(), "WeinR_intermediate")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (fpath in qc_obj@files) {
    
    key <- basename(fpath)
    if (verbose) message("\n-> Processing file: ", fpath)
    
    single_qc <- qc_obj
    single_qc@files <- fpath
    single_qc@metrics[[key]] <- qc_obj@metrics[[key]]
    single_qc@summary_metrics <- qc_obj@summary_metrics[
      qc_obj@summary_metrics$file == key, , drop = FALSE
    ]
    
    tmp1 <- NULL
    tmp2 <- NULL
    current_input <- fpath
    
    will_run_adapter <- !is.null(AdapterSeq)
    will_run_trim    <- !is.null(Start) || !is.null(End)
    
    failed <- FALSE
    skip_file <- FALSE
    
    single_qc <- tryCatch({
      
      # 1. Filter
      if (isTRUE(filter)) {
        if (verbose) message("   - Filtering...")
        
        is_final <- !(will_run_adapter || will_run_trim)
        
        single_qc <- FilterLong(
          qc_obj = single_qc,
          MinAvgQS = MinAvgQS,
          MinLength = MinLength,
          MaxNumberNs = MaxNumberNs,
          OutFileType = OutFileType,
          OutDir = processed_dir,
          WriteIntermediate = TRUE,
          KeepIntermediates = KeepIntermediates || is_final,
          OutSuffix = if (is_final) FinalSuffix else "filtered"
        )
        
        filt_sum <- single_qc@metadata[[key]]$filter_summary
        if (!is.null(filt_sum) && isTRUE(filt_sum$reads_after == 0L)) {
          if (verbose) message("   - Filtering kept 0 reads for ", key, ". Skipping remaining steps.")
          qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
          skip_file <- TRUE
        }
        
        if (!skip_file) {
          tmp1 <- attr(single_qc, "._tmp_fastq")
          if (!is.null(tmp1) && file.exists(tmp1)) {
            current_input <- tmp1
          } else if (!is_final) {
            if (verbose) message("   - No intermediate FASTQ from filtering for ", key, ". Skipping remaining steps.")
            qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
            skip_file <- TRUE
          }
        }
      }
      
      # 2. Adapter Trim
      if (!skip_file && !is.null(AdapterSeq)) {
        if (verbose) message("   - Adapter trimming...")
        
        is_final <- !will_run_trim
        
        single_qc <- RemoveAdapter(
          qc_obj = single_qc,
          adapterSeq = AdapterSeq,
          MaxMismatchEnd = MaxMismatchEnd,
          MinOverlapEnd = MinOverlapEnd,
          MinInternalDistance = MinInternalDistance,
          MinFragmentLength = MinFragmentLength,
          FilePath = current_input,
          OutDir = processed_dir,
          OutFileType = OutFileType,
          verbose = verbose,
          WriteIntermediate = TRUE,
          KeepIntermediates = KeepIntermediates || is_final,
          OutSuffix = if (is_final) FinalSuffix else "adaptertrimmed"
        )
        
        adap_sum <- single_qc@metadata[[key]]$adapter_summary
        if (!is.null(adap_sum) && isTRUE(adap_sum$reads_after == 0L)) {
          if (verbose) message("   - Adapter trimming kept 0 reads for ", key, ". Skipping remaining steps.")
          qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
          skip_file <- TRUE
        }
        
        if (!skip_file) {
          tmp2 <- attr(single_qc, "._tmp_fastq")
          if (!is.null(tmp2) && file.exists(tmp2)) {
            current_input <- tmp2
          }
        }
      }
      
      # 3. Start/End Trim
      if (!skip_file && (!is.null(Start) || !is.null(End))) {
        if (verbose) message("   - Start/end trimming...")
        
        single_qc <- TrimLong(
          qc_obj = single_qc,
          Start = Start,
          End = End,
          FilePath = current_input,
          OutDir = processed_dir,
          OutFileType = OutFileType,
          OutSuffix = FinalSuffix
        )
        
        trim_sum <- single_qc@metadata[[key]]$trim_summary
        if (!is.null(trim_sum) && isTRUE(trim_sum$reads_after == 0L)) {
          if (verbose) message("   - End trimming kept 0 reads for ", key, ".")
          qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
          skip_file <- TRUE
        }
      }
      
      single_qc
      
    }, error = function(e) {
      message("Failed [", key, "] - ", conditionMessage(e))
      failed <<- TRUE
      NULL
    })
    
    # If file skipped or failed, move on to next file
    if (failed || isTRUE(skip_file) || is.null(single_qc)) next
    
    qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
    
    # Optional cleanup of intermediates (if you later decide to enforce it here)
    if (!isTRUE(KeepIntermediates)) {
      if (!is.null(tmp1) && file.exists(tmp1)) unlink(tmp1)
      if (!is.null(tmp2) && file.exists(tmp2)) unlink(tmp2)
    }
  }
  
  message("\nProcessing complete for all files.")
  
  # 4. Report
  reports_dir <- file.path(base_dir, "Reports")
  dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(reports_dir, outpath)
  
  rmd_path  <- paste0(report_path, ".Rmd")
  html_path <- paste0(report_path, ".html")
  
  if (!force && (file.exists(rmd_path) || file.exists(html_path))) {
    stop(
      "Report already exists at:\n",
      rmd_path, "\n",
      "Use force = TRUE to overwrite.",
      call. = FALSE
    )
  }
  
  if (isTRUE(render_report)) {
    Reporter(
      mfa         = qc_obj,
      path        = report_path,
      title       = title,
      overwrite   = TRUE,
      render_html = TRUE,
      metadata    = TRUE
    )
  }
  
  qc_obj
}

