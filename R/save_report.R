

save_report <- function(qc_obj, file, base_dir, force, report_name) {
  
  print("Saving report...")
  report_path <- file.path(paste0(base_dir, "/reports/", report_name))
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
  Reporter(
    mfa         = qc_obj,
    path        = report_path,
    title       = report_name,
    overwrite   = TRUE,
    render_html = TRUE
  )
    
}