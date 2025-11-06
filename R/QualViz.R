#'QualViz is a one stop shop facilitator function 
#'


QualViz <- function(files, outpath, title){
  
  qc_obj <-.init_qc_object(files)
  
  for (file in qc_obj@files) {
    qsds <- ImportFile(file)[[1]]      # <-- pull the QSDS out of the list
    qc_obj <- QualMat(qc_obj, qsds, basename(file))
  }
  # (optional) make sure the output dir exists + use absolute path
  outpath <- normalizePath(outpath, winslash = "/", mustWork = FALSE)
  dir.create(dirname(outpath), recursive = TRUE, showWarnings = FALSE)
  
  CreateReport(mfa = qc_obj, path = outpath, title = title,
               overwrite = TRUE, render_html = TRUE)
  
  qc_obj
  
  #Return the LongReadQC object
  #base::return(qc_obj)
  
}