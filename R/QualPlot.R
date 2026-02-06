#' Generate and save QC plots for a single file
#'
#' Builds QC plots for one `filename` found in `qc_obj@metrics`, writes each plot
#' to disk as a high-resolution PNG, and stores the resulting file paths in
#' `qc_obj@plots[[filename]]`.
#'
#' @param qc_obj A `LongReadQC` object with populated `@metrics`.
#' @param filename Character(1). Must match a name in `names(qc_obj@metrics)`.
#' @param out_dir Character(1). Directory where plot outputs are written.
#'   A per-file subdirectory is created under `out_dir` using `basename(filename)`
#'   (sanitized).
#' @param width Numeric(1). Plot width in inches.
#' @param height Numeric(1). Plot height in inches.
#' @param dpi Integer(1). PNG resolution (DPI).
#' @param overwrite Logical(1). If TRUE, overwrite existing plot files.
#'
#' @return Updated `qc_obj` with plot file paths stored in `qc_obj@plots[[filename]]`.
#' @keywords internal
QualPlot <- function(qc_obj,
                     filename,
                     out_dir = "plots",
                     width = 7,
                     height = 5,
                     dpi = 300,
                     overwrite = TRUE) {
  
  # ---- Validate filename ----
  if (!base::is.character(filename) || base::length(filename) != 1L) {
    base::stop("`filename` must be a length-1 character string.")
  }
  
  if (base::is.null(qc_obj@metrics) || !base::is.list(qc_obj@metrics)) {
    base::stop("`qc_obj@metrics` must be a non-NULL list.")
  }
  
  if (!filename %in% base::names(qc_obj@metrics)) {
    base::stop(base::sprintf("File '%s' not found in qc_obj@metrics.", filename))
  }
  
  metrics <- qc_obj@metrics[[filename]]
  if (base::is.null(metrics) || !base::is.list(metrics)) {
    base::stop(base::sprintf("`qc_obj@metrics[['%s']]` must be a non-NULL list.", filename))
  }
  
  needed <- c("readLengths", "meanprQscore", "prGCcontent", "perPosQuality")
  missing <- needed[!needed %in% base::names(metrics)]
  if (base::length(missing) > 0L) {
    base::stop(base::sprintf(
      "Missing metrics for '%s': %s",
      filename,
      base::paste(missing, collapse = ", ")
    ))
  }
  
  # ---- Initialize plots slot if needed ----
  if (base::is.null(qc_obj@plots) || !base::is.list(qc_obj@plots)) {
    qc_obj@plots <- list()
  }
  
  # ---- Output directories ----
  base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # per-file folder name (sanitized, no extension)
  file_key <- base::gsub(
    "[^A-Za-z0-9._-]+", "_",
    tools::file_path_sans_ext(base::basename(filename))
  )
  
  file_dir <- base::file.path(out_dir, file_key)
  base::dir.create(file_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Helper: save PNG then discard plot ----
  save_png <- function(p, stem) {
    if (base::is.null(p)) return(NA_character_)
    
    out_path <- base::file.path(file_dir, base::paste0(stem, ".png"))
    if (base::file.exists(out_path) && !overwrite) return(out_path)
    
    grDevices::png(
      filename = out_path,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      type = "cairo-png",
      bg = "white"
    )
    
    # ensure device closes even if printing fails
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    
    print(p)
    
    # device close happens via on.exit
    out_path
  }
  
  # ---- Build → save → discard ----
  out_paths <- list(
    length_hist  = save_png(make_length_plot(metrics$readLengths), "length_hist"),
    quality_hist = save_png(make_quality_plot(metrics$meanprQscore), "quality_hist"),
    gc_hist      = save_png(make_gc_plot(metrics$prGCcontent), "gc_hist"),
    q_vs_length  = save_png(
      make_quality_vs_length_plot(metrics$readLengths, metrics$meanprQscore),
      "q_vs_length"
    ),
    per_pos_q    = save_png(make_per_position_plot(metrics$perPosQuality), "per_pos_q")
  )
  
  qc_obj@plots[[filename]] <- out_paths
  
  base::message(base::sprintf(
    "Saved QC plots for %s -> %s",
    base::basename(filename),
    file_dir
  ))
  
  qc_obj
}
