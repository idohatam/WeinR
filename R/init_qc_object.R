#' This is an internal function that is called by the main visualization workflow function.
#' It initiates the LongReadQC objects and fills in the files and metadata slots. 
#' @keywords internal
#' @noRd
#' 

.init_qc_object <- function(files) {
  new("LongReadQC",
      files = files,
      metadata = list(
        analysis_date = Sys.time(),
        n_files = length(files)
      ))
}