#' Run QC on input files and (optionally) generate a report
#'
#' Runs quality control on one or more long-read FASTQ files, producing summary
#' metrics and plots stored in a `LongReadQC` object. Optionally generates an
#' HTML report and saves the QC object as an `.rds` file.
#'
#' Outputs are written to:
#' \itemize{
#'   \item \code{WeinR_Outputs/Reports/} for \code{.Rmd} and \code{.html}
#'   \item \code{WeinR_Outputs/qc.rds} for the saved QC object
#' }
#'
#' @param files Character vector of input file paths.
#' @param report_name Character. Base filename (without extension) for the report.
#' @param render_report Logical. If \code{TRUE}, generate the HTML report.
#' @param force Logical. If \code{TRUE}, overwrite existing report files.
#'
#' @return A `LongReadQC` object containing metrics, plots, and summary tables for
#'   each input file (also saved to \code{WeinR_Outputs/qc.rds}).
#'
#' @examples
#' \dontrun{
#' fastq_files <- c(
#'   "data/Nakazawaea_holstii.fastq",
#'   "data/Nakazawaea_populi.fastq"
#' )
#'
#' qc_obj <- LinkAndReport(
#'   files = fastq_files,
#'   report_name = "Fully_linked_test",
#'   render_report = TRUE,
#'   force = TRUE
#' )
#'
#' qc_obj2 <- readRDS(file.path(getwd(), "WeinR_Outputs", "qc.rds"))
#' }
#'
#' @export
LinkAndReport <- function(
    files,
    report_name,
    render_report = TRUE,
    force = FALSE
) {
  
  qc_obj <- .init_qc_object(files)
  
  WeinRdir <- file.path(getwd(), "WeinR_Outputs")
  dir.create(WeinRdir, recursive = TRUE, showWarnings = FALSE)
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  base_dir <- file.path(WeinRdir, paste0("run_", timestamp))
  
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  
  metrics_dir <- file.path(base_dir, "metrics")
  dir.create(metrics_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (file in qc_obj@files) {
    qsds  <- ImportFile(file)
    fname <- basename(file)
    qc_obj <- QualMat(qc_obj, qsds, fname)
    #qc_obj <- QualPlot(qc_obj, filename = fname)
    file_metrics_dir <- file.path(metrics_dir, fname)
    
    dir.create(file_metrics_dir, recursive = TRUE, showWarnings = FALSE)
    
    metrics_list <- qc_obj@metrics[[fname]]
    
    for (metric_name in names(metrics_list)) {
      out_csv <- file.path(file_metrics_dir, paste0(metric_name, ".csv"))
      
      if (!force && file.exists(out_csv)) {
        stop(
          "Metric CSV already exists at:\n",
          out_csv, "\n",
          "Use force = TRUE to overwrite.",
          call. = FALSE
        )
      }
      
      metric <- metrics_list[[metric_name]]
      
      # normalize to data.frame for write.csv
      if (is.data.frame(metric)) {
        df <- metric
      } else if (is.matrix(metric)) {
        df <- as.data.frame(metric)
      } else if (is.list(metric)) {
        df <- data.frame(value = unlist(metric))
      } else {
        df <- data.frame(value = metric)
      }
      
      utils::write.csv(df, file = out_csv, row.names = FALSE)
    }
  }
  
  sm <- qc_obj@summary_metrics
  
  if (!"file" %in% colnames(sm)) {
    stop("qc_obj@summary_metrics must contain a 'file' column.", call. = FALSE)
  }
  
  out_csv <- file.path(metrics_dir, paste0(report_name, "_summary_metrics.csv"))
  
  if (!force && file.exists(out_csv)) {
    stop(
      "Metrics CSV already exists at:\n",
      out_csv, "\n",
      "Use force = TRUE to overwrite.",
      call. = FALSE
    )
  }

  utils::write.csv(sm, file = out_csv, row.names = FALSE)
  
  reports_dir <- file.path(base_dir, "Reports")
  dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
  
  report_path <- file.path(reports_dir, report_name)
  
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
    CreateReport(
      mfa         = qc_obj,
      path        = report_path,
      title       = report_name,
      overwrite   = TRUE,
      render_html = TRUE
    )
  }
  
  rds_path <- file.path(base_dir, "qc.rds")
  saveRDS(qc_obj, rds_path)
  
  qc_obj
}
