#' Read quality vs read length (rasterized scatter)
#'
#' Produces a scatter plot of per-read mean Phred Q-score vs read length,
#' with a log10 x-axis. Points are rasterized via ggrastr for performance.
#'
#' @param readLengths Numeric vector/list of read lengths.
#' @param meanprQscore Numeric vector/list of per-read mean Phred Q-scores.
#'
#' @return A ggplot object, or NULL if no valid points are available.
#' @keywords internal
make_quality_vs_length_plot <- function(readLengths, meanprQscore, file_name) {
  
  df <- data.frame(
    read_length = as.numeric(base::unlist(readLengths, use.names = FALSE)),
    avg_q_score = as.numeric(base::unlist(meanprQscore, use.names = FALSE))
  )
  
  # Clean for log10 scale + avoid warnings
  keep <- base::is.finite(df$read_length) &
    base::is.finite(df$avg_q_score) &
    (df$read_length > 0)
  
  df <- df[keep, , drop = FALSE]
  
  if (base::nrow(df) == 0L) {
    return(NULL)
  }
  
  ggplot2::ggplot(df, ggplot2::aes(x = read_length, y = avg_q_score)) +
    ggrastr::geom_point_rast(
      alpha = 0.25,
      color = "#0072B2",
      size = 0.4
    ) +
    
    # LOG scale fixes compression of dense region
    ggplot2::scale_x_log10(
      labels = scales::comma_format(),
      expand = c(0, 0)
    ) +
    
    ggplot2::scale_y_continuous(
      limits = c(0, 45),
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 10)
    ) +
    
    ggplot2::labs(
      x = "Read Length (bp, log scale)",
      y = "Average Q-Score",
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
