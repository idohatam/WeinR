#' Per-position quality plot
#'
#' Bins per-position quality values by read position, summarizes quantiles per bin,
#' applies loess smoothing to each quantile curve, and plots ribbons with a median line.
#'
#' @param perPosQuality Data.frame with columns `position` and `mean`.
#' @param window Integer bin size in bp.
#'
#' @return A ggplot object, or NULL if input is empty.
#' @keywords internal
make_per_position_plot <- function(perPosQuality, window = 2000) {
  
  # Input checks
  if (base::is.null(perPosQuality) || base::nrow(perPosQuality) == 0L) {
    return(NULL)
  }
  
  if (!("position" %in% base::names(perPosQuality)) ||
      !("mean" %in% base::names(perPosQuality))) {
    stop("`perPosQuality` must contain columns named 'position' and 'mean'.")
  }
  
  # Bin positions
  perPosQuality$bin <- base::floor(perPosQuality$position / window)
  
  # Summarise per bin
  df <- perPosQuality |>
    dplyr::group_by(.data$bin) |>
    dplyr::summarise(
      q10    = stats::quantile(.data$mean, 0.10, na.rm = TRUE),
      q25    = stats::quantile(.data$mean, 0.25, na.rm = TRUE),
      median = stats::median(.data$mean, na.rm = TRUE),
      q75    = stats::quantile(.data$mean, 0.75, na.rm = TRUE),
      q90    = stats::quantile(.data$mean, 0.90, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (base::nrow(df) == 0L) {
    return(NULL)
  }
  
  df$pos <- df$bin * window
  
  # ---- Inline loess smoothing ----
  fit_median <- stats::loess(median ~ pos, data = df, span = 0.3)
  fit_q25    <- stats::loess(q25    ~ pos, data = df, span = 0.3)
  fit_q75    <- stats::loess(q75    ~ pos, data = df, span = 0.3)
  fit_q10    <- stats::loess(q10    ~ pos, data = df, span = 0.3)
  fit_q90    <- stats::loess(q90    ~ pos, data = df, span = 0.3)
  
  df$median_s <- stats::predict(fit_median)
  df$q25_s    <- stats::predict(fit_q25)
  df$q75_s    <- stats::predict(fit_q75)
  df$q10_s    <- stats::predict(fit_q10)
  df$q90_s    <- stats::predict(fit_q90)
  
  # ---- Plot ----
  ggplot2::ggplot(df, ggplot2::aes(x = pos)) +
    
    # Outer 10–90% band
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = q10_s, ymax = q90_s),
      fill = "grey70",
      alpha = 0.4,
      color = NA
    ) +
    
    # Inner 25–75% band
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = q25_s, ymax = q75_s),
      fill = "#56B4E9",
      alpha = 0.6,
      color = NA
    ) +
    
    # Median line
    ggplot2::geom_line(
      ggplot2::aes(y = median_s),
      color = "#0072B2",
      linewidth = 1.4
    ) +
    
    ggplot2::labs(
      # title = "Per-Position Quality",
      subtitle = base::sprintf("Binned at %d bp", window),
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
      limits = c(0, base::max(df$pos, na.rm = TRUE) + 10000),
      labels = function(x) base::paste0(x / 1000, "k"),
      breaks = scales::extended_breaks(n = 15),
      expand = c(0, 0)
    ) +
    
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 10),
      expand = c(0, 0)
    )
}
