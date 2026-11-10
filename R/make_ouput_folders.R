

make_output_folders <- function(qc_obj) {
  
  WeinRdir <- file.path(getwd(), "WeinR_Outputs")
  dir.create(WeinRdir, recursive = TRUE, showWarnings = FALSE)
  
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  base_dir <- file.path(WeinRdir, paste0("run_", timestamp))
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  
  metrics_dir <- file.path(base_dir, "metrics")
  dir.create(metrics_dir, recursive = TRUE, showWarnings = FALSE)
  
  plots_dir <- file.path(base_dir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(file in qc_obj@files) {
    filename <- basename(file)
    file_metrics_dir <- file.path(metrics_dir, filename)
    dir.create(file_metrics_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  objects_dir <- file.path(base_dir, "objects")
  dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
  
  reports_dir <- file.path(base_dir, "reports")
  dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
  
  base_dir
}