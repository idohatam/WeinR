#' Process long-read files through filtering, adapter removal, end trimming, and reporting
#'
#' Runs a single end-to-end workflow over all files in a `LongReadQC` object. For each
#' file, this function can (optionally) apply read filtering (`FilterLong()`), adapter
#' trimming (`RemoveAdapter()`), and 5'/3' end trimming (`TrimLong()`), chaining steps
#' via temporary FASTQ files. After processing, it can generate an HTML report via
#' `CreateReport()`.
#'
#' @param qc_obj A `LongReadQC` object containing one or more input files in `qc_obj@files`,
#'   with corresponding metrics already computed in `qc_obj@metrics` and
#'   `qc_obj@summary_metrics`.
#' @param filter Logical(1). If `TRUE`, run `FilterLong()` prior to other steps.
#' @param MinAvgQS Numeric(1). Minimum mean per-read quality score to keep a read (passed to `FilterLong()`).
#' @param MinLength Integer(1). Minimum read length to keep a read (passed to `FilterLong()`).
#' @param MaxNumberNs Integer(1). Maximum number of `N` bases allowed per read (passed to `FilterLong()`).
#' @param AdapterSeq Character(1) or `NULL`. If not `NULL`, run `RemoveAdapter()` using this adapter sequence.
#' @param MaxMismatchEnd Integer(1). Maximum mismatches allowed in end-adapter matches (passed to `RemoveAdapter()`).
#' @param MinOverlapEnd Integer(1). Minimum overlap length for end-adapter matches (passed to `RemoveAdapter()`).
#' @param MinInternalDistance Integer(1). Minimum distance between internal adapter sites (passed to `RemoveAdapter()`).
#' @param MinFragmentLength Integer(1). Minimum fragment length to keep after internal adapter trimming (passed to `RemoveAdapter()`).
#' @param Start Integer(1) or `NULL`. If not `NULL`, trim this many bases from the 5' end (passed to `TrimLong()`).
#' @param End Integer(1) or `NULL`. If not `NULL`, trim this many bases from the 3' end (passed to `TrimLong()`).
#' @param OutFileType Character(1). Output file type for written reads (currently `"fastq"`).
#' @param outpath Character(1). Basename for the report output (without extension). The report is written to
#'   `file.path(getwd(), "reports", outpath)`, producing `.Rmd` and `.html`.
#' @param title Character(1). Title used in the generated report. Defaults to `outpath`.
#' @param render_report Logical(1). If `TRUE`, generate an HTML report with `CreateReport()`.
#' @param force Logical(1). If `FALSE` and an existing report `.Rmd` or `.html` already exists for `outpath`,
#'   the function errors. If `TRUE`, allows overwrite.
#' @param verbose Logical(1). If `TRUE`, print progress messages for each file and step.
#'
#' @return The updated `LongReadQC` object with per-file metadata updated in `qc_obj@metadata`.
#'   Processed reads are written to `OutDir` by the underlying step functions.
#'
#' @details
#' **Workflow order (per file):**
#' \enumerate{
#'   \item Optional filtering with `FilterLong()` if `filter = TRUE`.
#'   \item Optional adapter trimming with `RemoveAdapter()` if `AdapterSeq` is provided.
#'   \item Optional 5'/3' trimming with `TrimLong()` if `Start` and/or `End` are provided.
#' }
#'
#' Steps are chained using temporary FASTQ files stored as `attr(single_qc, "._tmp_fastq")`
#' produced by `FilterLong()` / `RemoveAdapter()` when `WriteIntermediate = TRUE`. Temporary
#' files are cleaned up automatically after each file.
#'
#' If filtering is enabled and no reads pass the filters, remaining steps for that file are
#' skipped and the workflow continues to the next input file.
#'
#' Reports are written to a `reports/` directory under the current working directory.
#' The report path is checked before rendering; set `force = TRUE` to overwrite.
#'
#' @seealso [FilterLong()], [RemoveAdapter()], [TrimLong()], [CreateReport()]
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
                         MinLength = 100,
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
                         verbose = TRUE) {
  
  message("Starting processing workflow for ", length(qc_obj@files), " file(s)...")
  
  base_dir <- file.path(getwd(), "WeinR_Outputs")
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  
  filtered_dir <- file.path(base_dir, "Processed_Files")
  filtered_dir <- normalizePath(filtered_dir, winslash = "/", mustWork = FALSE)
  dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (fpath in qc_obj@files) {
    
    key <- basename(fpath)  # <-- IMPORTANT: basename key used by QualMat/LinkAndReport
    
    tryCatch({
      
      if (verbose) message("\n-> Processing file: ", fpath)
      
      # Build single-file subobject
      single_qc <- qc_obj
      single_qc@files <- fpath
      
      # Copy metrics/summary using basename key
      single_qc@metrics[[key]] <- qc_obj@metrics[[key]]
      single_qc@summary_metrics <- qc_obj@summary_metrics[
        qc_obj@summary_metrics$file == key, , drop = FALSE
      ]
      
      tmp1 <- NULL  # filter temp
      tmp2 <- NULL  # adapter temp
      
      # Always cleanup temps created for this file
      on.exit({
        if (!is.null(tmp1) && file.exists(tmp1)) unlink(tmp1)
        if (!is.null(tmp2) && file.exists(tmp2)) unlink(tmp2)
      }, add = TRUE)
      
      current_input <- NULL
      
      # ---- Filter ----
      if (isTRUE(filter)) {
        if (verbose) message("   - Filtering...")
        
        single_qc <- FilterLong(
          qc_obj = single_qc,
          MinAvgQS = MinAvgQS,
          MinLength = MinLength,
          MaxNumberNs = MaxNumberNs,
          OutFileType = OutFileType,
          OutDir = filtered_dir,
          WriteIntermediate = TRUE
        )
        
        tmp1 <- attr(single_qc, "._tmp_fastq")
        
        # nothing passed filters -> skip remaining steps for this file
        if (is.null(tmp1) || !file.exists(tmp1)) {
          message("No reads passed filters for: ", key, ". Skipping remaining steps.")
          qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
          return(NULL)
        }
        
        current_input <- tmp1
      }
      
      # ---- Adapter Trim ----
      if (!is.null(AdapterSeq)) {
        if (verbose) message("   - Adapter trimming...")
        
        single_qc <- RemoveAdapter(
          qc_obj = single_qc,
          adapterSeq = AdapterSeq,
          MaxMismatchEnd = MaxMismatchEnd,
          MinOverlapEnd = MinOverlapEnd,
          MinInternalDistance = MinInternalDistance,
          MinFragmentLength = MinFragmentLength,
          FilePath = current_input,
          OutDir = filtered_dir,
          OutFileType = OutFileType,
          OutFile = NULL,
          verbose = verbose,
          WriteIntermediate = TRUE
        )
        
        tmp2 <- attr(single_qc, "._tmp_fastq")
        
        if (!is.null(tmp2) && file.exists(tmp2)) {
          current_input <- tmp2
        } else {
          current_input <- NULL
        }
      }
      
      # ---- 5'/3' end trim ----
      if (!is.null(Start) || !is.null(End)) {
        if (verbose) message("   - Start/end trimming...")
        
        single_qc <- TrimLong(
          qc_obj = single_qc,
          Start = Start,
          End = End,
          FilePath = current_input,
          OutDir = filtered_dir,
          OutFile = NULL,
          OutFileType = OutFileType
        )
      }
      
      # Store metadata back in main object (basename key)
      qc_obj@metadata[[key]] <- single_qc@metadata[[key]]
      
      NULL
      
    }, error = function(e) {
      
      warning(
        "Processing failed for file: ", key, "\n",
        "Reason: ", conditionMessage(e)
      )
      
      NULL
    })
  }
  
  message("\nProcessing complete for all files.")
  
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
      metadata = TRUE
    )
  }
  
  qc_obj
}

