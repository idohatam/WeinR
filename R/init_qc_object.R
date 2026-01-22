#' This is an internal function that is called by the main visualization workflow function.
#' It initiates the LongReadQC objects and fills in the files and metadata slots. 
#' @keywords internal
#' @noRd
#' 

.init_qc_object <- function(files) {
  new("LongReadQC",
      files = files,
      metadata = list(
        #analysis_date = Sys.time(),
        #n_files = length(files)
      ),
      summary_metrics = data.frame(
        file = character(length(files)),
        yield = numeric(length(files)),
        N50 = numeric(length(files)),
        N90 = numeric(length(files)),
        avgQscore = numeric(length(files)),
        `N count` = numeric(length(files)),
        stringsAsFactors = FALSE
      ))
}