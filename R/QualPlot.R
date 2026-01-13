#' Generate QC plots for a single file
#'
#' Builds and stores plots for one `filename` found in `qc_obj@metrics`.
#'
#' @param qc_obj An S4 LongReadQC object with `@metrics` and `@plots` slots.
#' @param filename Length-1 character string; must be a name in `qc_obj@metrics`.
#'
#' @return The updated `qc_obj` with plots stored in `qc_obj@plots[[filename]]`.
#' @keywords internal
#' 

QualPlot <- function(qc_obj, filename) {
  
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
  
  # ---- Optional: ensure expected metric components exist ----
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
  
  # ---- Build plots ----
  qc_obj@plots[[filename]] <- list(
    length_hist  = make_length_plot(metrics$readLengths),
    quality_hist = make_quality_plot(metrics$meanprQscore),
    gc_hist      = make_gc_plot(metrics$prGCcontent),
    q_vs_length  = make_quality_vs_length_plot(
      metrics$readLengths,
      metrics$meanprQscore
    ),
    per_pos_q    = make_per_position_plot(metrics$perPosQuality)
  )
  
  base::message(base::sprintf("Processed file: %s", filename))
  qc_obj
}

# ==== How to use: Example Usage ====

#qc_obj <- QualPlot(qc_obj, filename = "sample1.fastq")

# Visualize
# print(qc_obj@plots[["sample1.fastq"]]$length_hist)
# print(qc_obj@plots[["sample1.fastq"]]$quality_hist)
# print(qc_obj@plots[["sample1.fastq"]]$gc_hist)
#print(qc_obj@plots[["sample1.fastq"]]$q_vs_length)
# print(qc_obj@plots[["sample1.fastq"]]$per_pos_q)
