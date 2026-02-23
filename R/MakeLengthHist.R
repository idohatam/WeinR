#' Read length histogram (log10 x-axis)
#'
#' Produces a histogram of read lengths with a log10 x-axis.
#'
#' @param readLengths Numeric vector (or list coercible to numeric) of read lengths.
#'
#' @return A ggplot object.
#' @keywords internal
make_length_plot <- function(readLengths, file_name) {
  
  readLengths <- as.numeric(unlist(readLengths, use.names = FALSE))
  readLengths <- readLengths[is.finite(readLengths) & readLengths > 0]
  
  ggplot2::ggplot(
    data.frame(read_length = readLengths),
    ggplot2::aes(x = read_length)
  ) +
    ggplot2::geom_histogram(
      bins = 50,
      fill = "#0072B2",
      color = "white",
      alpha = 0.9
    ) +
    ggplot2::scale_x_log10(
      labels = scales::comma_format(),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 8)
    ) +
    ggplot2::labs(
      x = "Read Length (bp, log scale)",
      y = "Read Count",
      title = file_name
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
      axis.line  = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.6)
    )
}
