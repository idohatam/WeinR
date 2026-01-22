#' Run QC on input files and (optionally) generate a report
#'
#' @param files Character vector of input file paths.
#' @param outpath Character. Base filename (without extension) for the report.
#' @param title Character. Report title. Defaults to `outpath`.
#' @param render_report Logical. If `TRUE`, generate the report.
#' @param force Logical. If `TRUE`, overwrite existing report files.
#'
#' @export
LinkAndReport <- function(
    files,
    outpath,
    title = outpath,
    render_report = TRUE,
    force = FALSE
) {
  
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
  
  ## ---- FORCE CHECK ---------------------------------------------------------
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
  ## -------------------------------------------------------------------------
  
  if (isTRUE(render_report)) {
    CreateReport(
      mfa         = qc_obj,
      path        = report_path,
      title       = title,
      overwrite   = TRUE,   # controlled by force check above
      render_html = TRUE
    )
  }
  
  qc_obj
}
