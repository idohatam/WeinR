library(ggplot2)
library(Rsamtools)

#' Process per-read PolyA tail lengths
#'
#' Extracts PolyA tail lengths from the BAM file.
#'
#' @param bam_file Path to the BAM file.
#'
#' @return Numeric vector of valid PolyA tail lengths.
#' @keywords internal
process_polyA_lengths <- function(bam_file) {
  
  bam_data <- scanBam(
    bam_file,
    param = ScanBamParam(
      tag = "pt"
    )
  )
  
  polyA_lengths <- as.numeric(bam_data[[1]]$tag$pt)
  polyA_lengths <- polyA_lengths[is.finite(polyA_lengths) & polyA_lengths >= 0]
  
  polyA_lengths
}


#' PolyA tail length histogram
#'
#' Produces a histogram of per-read PolyA tail lengths.
#'
#' @param polyA_lengths Numeric vector of PolyA tail lengths.
#' @param file_name Character title for the plot.
#'
#' @return A ggplot object, or NULL if no valid PolyA lengths are available.
#' @keywords internal
make_polyA_plot <- function(polyA_lengths, file_name) {
  
  polyA_lengths <- as.numeric(unlist(polyA_lengths, use.names = FALSE))
  polyA_lengths <- polyA_lengths[is.finite(polyA_lengths) & polyA_lengths >= 0]
  
  if (length(polyA_lengths) == 0L) {
    return(NULL)
  }
  
  ggplot(
    data.frame(PolyA_length = polyA_lengths),
    aes(x = PolyA_length)
  ) +
    geom_histogram(
      bins = 50,
      fill = "#0072B2",
      color = "white",
      alpha = 0.9
    ) +
    scale_x_continuous(
      labels = scales::comma_format(),
      breaks = scales::pretty_breaks(n = 12),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 12)
    ) +
    labs(
      x = "PolyA Tail Length (nt)",
      y = "Number of Reads",
      title = file_name
    ) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6)
    )
}


# #bam_file <- "/Users/aaronsantos/Downloads/subsampled_1percent.bam"
# 
# polyA_lengths <- process_polyA_lengths(bam_file)
# 
# print(
#   make_polyA_plot(
#     polyA_lengths,
#     "PolyA Tail Length Distribution"
#   )
# )