# File-by-File Mode: QualPlot() Definition
# Compatible with LongReadQC structure

# Required libraries
#library(ggplot2)
#library(scales)
# library(mgcv) # for geom_smooth(method = "gam")
#library(ggrastr)
#pip install 
#

#break up into own file. 
# qual plot calls out helpers
# all other functions as own things
# work flow fuctions, so plotting functions arent helper function, 
# instead of 
# We already QualViz, so technically we have the one that runs thorugh and calls the inputs 
# Nicer alls helper functions as its own function, each helper gets its own
# don't need to state the packages
# try hex or geom_point_density and see if it works better than point_rast # ugly but works .
# all be in a plot folder


# Main function: process one file at a time
QualPlot <- function(qc_obj, filename) {
  
  if (!filename %in% names(qc_obj@metrics)) {
    stop(paste("File", filename, "not found in qc_obj@metrics"))
  }
  
  metrics <- qc_obj@metrics[[filename]]
  
  # Timing wrappers
  time_it <- function(label, expr) {
    start <- Sys.time()
    result <- eval(expr)
    end <- Sys.time()
    message(sprintf("%s for %s took %.3f seconds",
                    label, filename, as.numeric(end - start)))
    return(result)
  }
  
  qc_obj@plots[[filename]] <- list(
    length_hist  = time_it("make_length_plot()", quote(make_length_plot(metrics$readLengths))),
    quality_hist = time_it("make_quality_plot()", quote(make_quality_plot(metrics$meanprQscore))),
    gc_hist      = time_it("make_gc_plot()", quote(make_gc_plot(metrics$prGCcontent))),
    q_vs_length  = time_it("make_quality_vs_length_plot()", 
                           quote(make_quality_vs_length_plot(metrics$readLengths, metrics$meanprQscore))),
    per_pos_q    = time_it("make_per_position_plot()", 
                           quote(make_per_position_plot(metrics$perPosQuality)))
  )
  
  
  message(sprintf("Processed file: %s", filename))
  return(qc_obj)
}


# ==== Helper Plot Functions ====

# 1. Read Length Histogram (log scale)
make_length_plot <- function(readLengths) {
  
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
      title = "Read Length Distribution"
    ) +
    
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
      axis.line     = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks    = ggplot2::element_line(color = "black", linewidth = 0.6)
      
      # plot.margin = ggplot2::margin(t = 5, r = 10, b = 20, l = 20)
      
    )
}



# 2. Average Q Score Histogram with vertical lines
make_quality_plot <- function(meanprQscore, avgQScore = NULL) {
  
  meanprQscore <- unlist(meanprQscore)
  
  # Compute statistics
  mean_q   <- if (!is.null(avgQScore)) avgQScore else base::mean(meanprQscore, na.rm = TRUE)
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
      # ggplot2::aes(xintercept = mean_q, color = "Mean"),   # <-- legend mapping removed
      xintercept = mean_q,                                   # <-- direct color
      color = "#BF616A",
      linetype = "dashed",
      linewidth = 1.5
    ) +
    
    # Median line
    ggplot2::geom_vline(
      # ggplot2::aes(xintercept = median_q, color = "Median"), # <-- legend mapping removed
      xintercept = median_q,
      color = "#0072B2",
      linetype = "solid",
      linewidth = 1.5
    ) +
    
    # Legend palette (DISABLED)
    # ggplot2::scale_color_manual(
    #   name = NULL,
    #   values = c(
    #     "Mean"   = "#BF616A",
    #     "Median" = "#0072B2"
    #   )
    # ) +
    
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
      
      # Legend settings (DISABLED)
      legend.position = "none",
      # legend.position.inside = c(0.02, 0.96),
      # legend.justification = c(0, 1),
      # legend.direction = "horizontal",
      # legend.text = ggplot2::element_text(size = 12),
      # legend.key = ggplot2::element_blank(),
      
      plot.background  = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      
      plot.margin = ggplot2::margin(12, 15, 15, 15)
    )
  
  # Legend override settings (DISABLED)
  # +
  # ggplot2::guides(
  #   color = ggplot2::guide_legend(
  #     override.aes = list(
  #       linetype = c("dashed", "solid"),
  #       linewidth = c(2, 3)
  #     )
  #   )
  # )
}





# 3. GC Content Histogram (average per read)
make_gc_plot <- function(prGCContent) {
  
  # Flatten list safely
  prGCContent <- unlist(prGCContent)
  
  # Skip if empty
  if (length(prGCContent) == 0) return(NULL)
  
  # Compute maximum histogram bin height for proper y-axis limits
  # hist_vals <- graphics::hist(prGCContent, plot = FALSE, breaks = 40)$counts
  # upper_y <- ceiling(max(hist_vals) / 10) * 10
  
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
      title = "GC Content Distribution"
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


# 5. Q Score vs. Read Length (Scatter + Smooth) #geom point density
make_quality_vs_length_plot <- function(readLengths, meanprQscore) {
  
  df <- data.frame(
    read_length = as.numeric(unlist(readLengths, use.names = FALSE)),
    avg_q_score = as.numeric(unlist(meanprQscore, use.names = FALSE))
  )
  
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
      title = "Read Quality vs. Read Length"
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
      panel.background = ggplot2::element_blank(),
      
    )
}




# 6. Per-Position Quality Plot 
make_per_position_plot <- function(perPosQuality, window = 2000) {
  
  if (is.null(perPosQuality) || nrow(perPosQuality) == 0)
    return(NULL)
  
  perPosQuality$bin <- floor(perPosQuality$position / window)
  
  df <- perPosQuality |>
    dplyr::group_by(bin) |>
    dplyr::summarise(
      q10    = stats::quantile(mean, 0.10, na.rm = TRUE),
      q25    = stats::quantile(mean, 0.25, na.rm = TRUE),
      median = stats::median(mean, na.rm = TRUE),
      q75    = stats::quantile(mean, 0.75, na.rm = TRUE),
      q90    = stats::quantile(mean, 0.90, na.rm = TRUE)
    )
  
  df$pos <- df$bin * window
  
  # Smooth curves
  smooth_loess <- function(x, y) {
    fit <- stats::loess(y ~ x, span = 0.3)
    stats::predict(fit)
  }
  
  df$median_s <- smooth_loess(df$pos, df$median)
  df$q25_s    <- smooth_loess(df$pos, df$q25)
  df$q75_s    <- smooth_loess(df$pos, df$q75)
  df$q10_s    <- smooth_loess(df$pos, df$q10)
  df$q90_s    <- smooth_loess(df$pos, df$q90)
  
  ggplot2::ggplot(df, ggplot2::aes(x = pos)) +
    
    # Outer 10–90% band
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = q10_s,
        ymax = q90_s
        # , fill = "10–90% Range"
      ),
      fill = "grey70", alpha = 0.4, color = NA
    ) +
    
    # Inner 25–75% band
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = q25_s,
        ymax = q75_s
        # , fill = "25–75% IQR"
      ),
      fill = "#56B4E9", alpha = 0.6, color = NA
    ) +
    
    # Median line
    ggplot2::geom_line(
      ggplot2::aes(
        y = median_s
        # , color = "Median (Smoothed)"
      ),
      color = "#0072B2", linewidth = 1.4
    ) +
    
    ggplot2::labs(
      title = "Per-Position Quality",
      subtitle = "Binned at 2000 bp",
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
      limits = c(0, max(df$pos) + 10000),
      labels = function(x) paste0(x / 1000, "k"),
      breaks = scales::extended_breaks(n = 15),
      expand = c(0, 0)
    ) +
    
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 10),
      expand = c(0, 0)
    )
}



# ==== How to use: Example Usage ====

#qc_obj <- QualPlot(qc_obj, filename = "sample1.fastq")

# Visualize
# print(qc_obj@plots[["sample1.fastq"]]$length_hist)
# print(qc_obj@plots[["sample1.fastq"]]$quality_hist)
# print(qc_obj@plots[["sample1.fastq"]]$gc_hist)
#print(qc_obj@plots[["sample1.fastq"]]$q_vs_length)
# print(qc_obj@plots[["sample1.fastq"]]$per_pos_q)
