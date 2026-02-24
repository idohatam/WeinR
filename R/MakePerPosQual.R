#' Per-position quality plot (no binning)
#'
#' Applies rolling-mean smoothing (±2500 positions) to each quantile curve
#' and plots ribbons with a median line.
#'
#' @param perPosQuality Data.frame with columns: position, mean, median, q25, q75
#' @param file_name Character title for the plot.
#'
#' @return A ggplot object, or NULL if input is empty.
#' @keywords internal
make_per_position_plot <- function(perPosQuality, file_name) {
  
  # ---- Input checks ----
  if (base::is.null(perPosQuality) || base::nrow(perPosQuality) == 0L) {
    return(NULL)
  }
  
  required_cols <- c("position", "mean", "median", "q25", "q75")
  if (!all(required_cols %in% base::names(perPosQuality))) {
    stop("`perPosQuality` must contain columns: position, mean, median, q25, q75.")
  }
  
  # ---- Order by position ----
  df <- perPosQuality[base::order(perPosQuality$position), , drop = FALSE]
  
  # ---- Rolling smoothing (±2500 positions) ----
  half_window <- 2500L
  k <- 2L * half_window + 1L  # 5001
  
  
  if (base::nrow(df) > k) {
    dt <- data.table::as.data.table(df)
    n <- base::nrow(dt)
    
    # compute smoothed vectors
    mean_s   <- data.table::frollmean(dt[["mean"]],   n = k, align = "center")
    median_s <- data.table::frollmean(dt[["median"]], n = k, align = "center")
    q25_s    <- data.table::frollmean(dt[["q25"]],    n = k, align = "center")
    q75_s    <- data.table::frollmean(dt[["q75"]],    n = k, align = "center")
    
    # fill edges with value at index 2501 and n-2500
    mean_s[1:half_window]   <- mean_s[half_window + 1L]
    median_s[1:half_window] <- median_s[half_window + 1L]
    q25_s[1:half_window]    <- q25_s[half_window + 1L]
    q75_s[1:half_window]    <- q75_s[half_window + 1L]
    
    mean_s[(n - half_window + 1L):n]   <- mean_s[n - half_window]
    median_s[(n - half_window + 1L):n] <- median_s[n - half_window]
    q25_s[(n - half_window + 1L):n]    <- q25_s[n - half_window]
    q75_s[(n - half_window + 1L):n]    <- q75_s[n - half_window]
    
    dt[["mean"]]   <- mean_s
    dt[["median"]] <- median_s
    dt[["q25"]]    <- q25_s
    dt[["q75"]]    <- q75_s
    
    df <- base::as.data.frame(dt)
  }
  
  # ---- Plot ----
  ggplot2::ggplot(df, ggplot2::aes(x = position)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = q25, ymax = q75),
      fill = "#56B4E9",
      alpha = 0.6,
      color = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = median),
      color = "#0072B2",
      linewidth = 1.4
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = mean),
      color = "#007",
      linewidth = 1.4
    ) +
    ggplot2::labs(
      title = file_name,
      x = "Read Position (bp)",
      y = "Phred Q-Score"
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.line          = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks         = ggplot2::element_line(color = "black", linewidth = 0.6),
      plot.title         = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle      = ggplot2::element_text(hjust = 0),
      legend.position    = "none"
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) base::paste0(x / 1000, "k"),
      breaks = scales::extended_breaks(n = 15),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 10),
      expand = c(0, 0)
    )
}