#' QualViz is a one-stop facilitator:
#' - builds the LongReadQC object
#' - imports reads
#' - computes metrics (QualMat)
#' - makes plots (QualPlot)
#' - optionally writes a report (CreateReport)
LinkAndReport <- function(files, outpath, title, render_report = TRUE) {
  # init object
  qc_obj <- .init_qc_object(files)
  
  for (file in qc_obj@files) {
    qsds <- ImportFile(file)[[1]]     
    
    fname <- basename(file)
    qc_obj <- QualMat(qc_obj, qsds, fname)
    
    qc_obj <- QualPlot(qc_obj, filename = fname)
  }
  
  # ensure output location exists
  outpath <- normalizePath(outpath, winslash = "/", mustWork = FALSE)
  dir.create(dirname(outpath), recursive = TRUE, showWarnings = FALSE)
  
  if (isTRUE(render_report)) {
    CreateReport(
      mfa         = qc_obj,
      path        = outpath,
      title       = title,
      overwrite   = TRUE,
      render_html = TRUE
    )
  }

  return(qc_obj)
}



