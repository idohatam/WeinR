#' GC content histogram (average per read)
#'
#' Produces a histogram of per-read GC content values.
#'
#' @param prGCContent Numeric vector/list of per-read GC content values.
#'
#' @return A ggplot object, or NULL if no finite GC values are available.
#' @keywords internal
make_gc_plot <- function(prGCContent, file_name) {
  
  # Flatten and clean
  prGCContent <- as.numeric(base::unlist(prGCContent, use.names = FALSE))
  prGCContent <- prGCContent[base::is.finite(prGCContent)]
  
  # Skip if empty
  if (base::length(prGCContent) == 0L) {
    return(NULL)
  }
  
  ggplot2::ggplot(
    data = data.frame(gc_content = prGCContent),
    ggplot2::aes(x = gc_content)
  ) +
    ggplot2::geom_histogram(
      bins = 40,
      fill = "#CC79A7",
      color = "white",
      alpha = 0.9
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 10)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 10)
    ) +
    ggplot2::labs(
      x = "GC Content",
      y = "Read Count",
      title = file_name
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      
      axis.line  = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.6),
      
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      
      plot.background  = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    )
}
