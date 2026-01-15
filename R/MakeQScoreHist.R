#' Average Q-score histogram (with mean/median lines)
#'
#' Produces a histogram of per-read mean Phred Q-scores with vertical lines for
#' mean (dashed) and median (solid).
#'
#' @param meanprQscore Numeric vector/list of per-read mean Phred Q-scores.
#' @param avgQScore Optional numeric scalar. If provided, used as the mean line
#'   instead of computing the mean from `meanprQscore`.
#'
#' @return A ggplot object, or NULL if no finite Q-scores are available.
#' @keywords internal
make_quality_plot <- function(meanprQscore, avgQScore = NULL) {
  
  meanprQscore <- as.numeric(base::unlist(meanprQscore, use.names = FALSE))
  meanprQscore <- meanprQscore[base::is.finite(meanprQscore)]
  
  if (base::length(meanprQscore) == 0L) {
    return(NULL)
  }
  
  # Compute statistics
  mean_q <- if (!base::is.null(avgQScore)) {
    as.numeric(avgQScore)[1L]
  } else {
    base::mean(meanprQscore, na.rm = TRUE)
  }
  
  median_q <- stats::median(meanprQscore, na.rm = TRUE)
  
  ggplot2::ggplot(
    data = data.frame(meanprQscore = meanprQscore),
    ggplot2::aes(x = meanprQscore)
  ) +
    
    # Histogram
    ggplot2::geom_histogram(
      bins = 40,
      fill = "#8FBCBB",
      color = "white",
      alpha = 0.8
    ) +
    
    # Mean line
    ggplot2::geom_vline(
      xintercept = mean_q,
      color = "#BF616A",
      linetype = "dashed",
      linewidth = 1.5
    ) +
    
    # Median line
    ggplot2::geom_vline(
      xintercept = median_q,
      color = "#0072B2",
      linetype = "solid",
      linewidth = 1.5
    ) +
    
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 10)
    ) +
    
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 8)
    ) +
    
    ggplot2::labs(
      x = "Average Phred Q-Score",
      y = "Read Count",
      title = "Average Q-Score Distribution"
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
      
      legend.position = "none",
      
      plot.background  = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      
      plot.margin = ggplot2::margin(12, 15, 15, 15)
    )
}
