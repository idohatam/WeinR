#' Run QC on input files and return a LongReadQC object
#'
#' Internal convenience wrapper that initializes a `LongReadQC` object,
#' imports each input file, computes metrics, and generates plots.
#' This function returns the populated `LongReadQC` object and does not
#' render a report (report rendering is handled elsewhere).
#'
#' @param files Character vector of input file paths (e.g., FASTQ files).
#' @param outpath Character(1). Output base path (used to create the output
#'   directory if needed). Currently not used to write files in this function.
#' @param title Character(1). Title string (currently unused in this function).
#'
#' @return A `LongReadQC` object containing computed metrics and plots.
#'
#' @details
#' This function is intended for internal use and is not part of the public API.
#' It assumes `files` are readable by `ImportFile()` and that the internal helpers
#' `.init_qc_object()`, `QualMat()`, and `QualPlot()` exist.
#'
#' @examples
#' \dontrun{
#' qc <- WeinR:::QualViz(files = fastq_files, outpath = "reports/qc", title = "QC")
#' }
#'
#' @keywords internal
QualViz <- function(files, outpath, title){
  
  qc_obj <-.init_qc_object(files)
  
  for (file in qc_obj@files) {
    reads <- ImportFile(file)      
    qc_obj <- QualMat(qc_obj, reads, basename(file))
    qc_obj <- QualPlot(qc_obj, basename(file))
  }
  # (optional) make sure the output dir exists + use absolute path
  outpath <- normalizePath(outpath, winslash = "/", mustWork = FALSE)
  dir.create(dirname(outpath), recursive = TRUE, showWarnings = FALSE)
  
  #CreateReport(mfa = qc_obj, path = outpath, title = title,
  #             overwrite = TRUE, render_html = TRUE)
  
  #qc_obj
  
  #Return the LongReadQC object
  base::return(qc_obj)
  
}
