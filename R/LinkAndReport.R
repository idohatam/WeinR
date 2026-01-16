#' Run QC on input files and (optionally) generate a report
#'
#' Initializes a `LongReadQC` object, imports each file, computes metrics, and
#' generates plots. If `render_report = TRUE`, it calls [CreateReport()] to write
#' an `.Rmd` (and optionally render an `.html`) report into a `reports/` folder
#' under the current working directory.
#'
#' @param files Character vector of input file paths (e.g., FASTQ files).
#' @param outpath Character. Base filename (without extension) for the report.
#'   The report is written to `file.path(getwd(), "reports", outpath)`, and
#'   [CreateReport()] will append `.Rmd` if needed.
#' @param title Character. Report title. Defaults to `outpath`.
#' @param render_report Logical. If `TRUE`, generate the report via
#'   [CreateReport()]. If `FALSE`, only returns the QC object.
#'
#' @return A `LongReadQC` object containing computed metrics and plots.
#'
#' @details
#' This function assumes:
#' - `files` are readable by `ImportFile()`
#' - `.init_qc_object()`, `QualMat()`, and `QualPlot()` exist and return updated
#'   `LongReadQC` objects
#' - The report directory is `getwd()/reports`
#'
#' @examples
#' \dontrun{
#' qc <- LinkAndReport(
#'   files   = fastq_files,
#'   outpath = "Fully_linked_test"
#' )
#' }
#'
#' @export
LinkAndReport <- function(files, outpath, title = outpath, render_report = TRUE) {
  
  qc_obj <- .init_qc_object(files)
  
  for (file in qc_obj@files) {
    qsds  <- ImportFile(file)
    fname <- basename(file)
    qc_obj <- QualMat(qc_obj, qsds, fname)
    qc_obj <- QualPlot(qc_obj, filename = fname)
  }
  
  reports_dir <- file.path(getwd(), "reports")
  dir.create(reports_dir, showWarnings = FALSE, recursive = TRUE)
  
  report_path <- file.path(reports_dir, outpath)
  
  if (isTRUE(render_report)) {
    CreateReport(
      mfa         = qc_obj,
      path        = report_path,
      title       = title,
      overwrite   = TRUE,
      render_html = TRUE
    )
  }
  
  qc_obj
}