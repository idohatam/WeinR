#' LongReadQC Class
#'
#' An S4 class to store quality control results for long-read sequencing data.
#'
#' @importFrom methods setClass setMethod show new
#'
#' @slot files Character vector of input file paths.
#' @slot metrics Named list of data.frames containing per-read quality metrics.
#' @slot plots Nested list structure: file -> plot_type -> ggplot object.
#' @slot summary_metrics data.frame with aggregated statistics across files.
#' @slot metadata List containing analysis parameters, timestamps, and related information.
#'
#' @export
setClass("LongReadQC",
         slots = list(
           files = "character",
           metrics = "list",
           plots = "list",
           summary_metrics = "data.frame",
           metadata = "list"
         ),
         prototype = list(
           files = character(0),
           metrics = list(),
           plots = list(),
           summary_metrics = data.frame(),
           metadata = list()
         )
)

#' Show method for LongReadQC objects
#'
#' Displays a concise summary of a \code{LongReadQC} object, including the
#' number of processed files and whether metrics and plots are available.
#'
#' @param object A \code{LongReadQC} object.
#'
#' @return Invisibly returns the input \code{LongReadQC} object.
#'
#' @export
setMethod("show", "LongReadQC", function(object) {
  cat("LongReadQC object\n")
  cat("----------------\n")
  cat("Files processed:", length(object@files), "\n")
  if (length(object@files) > 0) {
    cat("  ", paste(basename(object@files), collapse = ", "), "\n")
  }
  cat("Metrics available:", length(object@metrics) > 0, "\n")
  cat("Plots available:", length(object@plots) > 0, "\n")
  invisible(object)
})