

save_results <- function(qc_obj, base_dir, force, report_name, render_report) {
  
  for(file in qc_obj@files) {
    qc_obj <- QualPlot(
      qc_obj,
      filename  = file,
      out_dir   = paste0(base_dir, "/plots"),
      dpi       = 300,
      overwrite = force
    )
    save_metrics(qc_obj, file, base_dir, force, report_name)
    save_object(qc_obj, file, base_dir)
    if(render_report) {
      save_report(qc_obj, file, base_dir, force, report_name)
    }
  }

}