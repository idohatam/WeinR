

save_object <- function(qc_obj, file, base_dir) {
  
  print("Saving object...")
  rds_path <- file.path(paste0(base_dir, "/objects"), "qc.rds")
  saveRDS(qc_obj, rds_path)
  
}