

save_metrics <- function(qc_obj, file, base_dir, force, report_name) {
  
  print("Saving metrics...")
  filename <- file
  metrics_list <- qc_obj@metrics[[filename]]
  for (metric_name in names(metrics_list)) {
    out_csv <- file.path(paste0(base_dir, "/metrics/", basename(filename)), paste0(metric_name, ".csv"))
    if (!force && file.exists(out_csv)) {
      stop(
        "Metric CSV already exists at:\n",
        out_csv, "\n",
        "Use force = TRUE to overwrite.",
        call. = FALSE
      )
    }
    metric <- metrics_list[[metric_name]]

    # normalize to data.frame for write.csv
    if (is.data.frame(metric)) {
      df <- metric
    } else if (is.matrix(metric)) {
      df <- as.data.frame(metric)
    } else if (is.list(metric)) {
      df <- data.frame(value = unlist(metric))
    } else {
      df <- data.frame(value = metric)
    }
    utils::write.csv(df, file = out_csv, row.names = FALSE)
    
    sm <- qc_obj@summary_metrics
  
    if (!"File" %in% colnames(sm)) {
      stop("qc_obj@summary_metrics must contain a 'File' column.", call. = FALSE)
    }
  
    out_csv <- file.path(paste0(base_dir, "/metrics"), paste0(report_name, "_summary_metrics.csv"))
  
    if (!force && file.exists(out_csv)) {
      stop(
        "Metrics CSV already exists at:\n",
        out_csv, "\n",
        "Use force = TRUE to overwrite.",
        call. = FALSE
      )
    }
  
    utils::write.csv(sm, file = out_csv, row.names = FALSE)
  }
  
}