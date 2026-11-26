library("Rcpp")
library("ShortRead")
library("stringr")
library("data.table")

create_report <- function(
    files,
    minimum_reads = 1L,
    report_name,
    render_report = TRUE,
    force = FALSE
) {
  
  qc_obj <- validate_input(files, minimum_reads)
  base_dir <- make_output_folders(qc_obj)
  for (file in qc_obj@files) {
    qc_obj <- process_file_report(file, qc_obj)
  }
  save_results(qc_obj, base_dir, force, report_name, render_report)

  qc_obj
}