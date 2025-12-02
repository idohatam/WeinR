# File-by-File Mode: QualPlot() Definition
# Compatible with LongReadQC structure

# Required libraries
# library(ggplot2)
# library(scales)
# library(mgcv) # for geom_smooth(method = "gam")
# library(ggrastr)pip install 
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
  ggplot(data.frame(read_length = readLengths), aes(x = read_length)) +
    geom_histogram(
      bins = 50,
      fill = "#0072B2",      
      color = "white",
      alpha = 0.9
    ) +
    scale_x_log10(labels = comma_format()) +
    labs(
      x = "Read Length (bp, log scale)",
      y = "Read Count",
      title = "Read Length Distribution"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),          # Removes all gridlines
      panel.background = element_blank(),    # Ensures clean background
      plot.background = element_blank(),     # Removes gray panel background
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.line = element_line(color = "black"),  # Adds clean x/y axis lines
      axis.ticks = element_line(color = "black")  # Adds small axis ticks
    )
}



# 2. Average Q Score Histogram with vertical lines
make_quality_plot <- function(meanprQscore, avgQScore = NULL) {
  meanprQscore <- unlist(meanprQscore)   # flatten list -> numeric vector
  
  mean_q <- if (!is.null(avgQScore)) avgQScore else mean(meanprQscore, na.rm = TRUE)
  median_q <- median(meanprQscore, na.rm = TRUE)
  
  ggplot(mapping = aes(x = meanprQscore)) +
    geom_histogram(
      bins = 40,
      fill = "#009E73",
      color = "white",
      alpha = 0.9
    ) +
    geom_vline(
      xintercept = mean_q,
      color = "#D55E00",
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_vline(
      xintercept = median_q,
      color = "#0072B2",
      linetype = "solid",
      linewidth = 1
    ) +
    labs(
      x = "Average Phred Q-Score",
      y = "Read Count",
      title = "Average Q-Score Distribution"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}



# 3. GC Content Histogram (average per read)
make_gc_plot <- function(prGCContent) {
  # Flatten list safely
  prGCContent <- unlist(prGCContent)
  
  # Skip if empty
  if (length(prGCContent) == 0) return(NULL)
  
  # Plot
  ggplot(data.frame(gc_content = prGCContent), aes(x = gc_content)) +
    geom_histogram(bins = 40, fill = "#CC79A7", color = "white", alpha = 0.9) +
    labs(
      x = "GC Content",
      y = "Read Count",
      title = "GC Content Distribution"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# 5. Q Score vs. Read Length (Scatter + Smooth) #geom point density
make_quality_vs_length_plot <- function(readLengths, meanprQscore) {
  df <- data.frame(
    read_length = as.numeric(unlist(readLengths, use.names = FALSE)),
    avg_q_score = as.numeric(unlist(meanprQscore, use.names = FALSE))
  )

  ggplot(df, aes(x = read_length, y = avg_q_score)) +
    geom_point_rast(alpha = 0.25, color = "#0072B2", size = 0.4) + # geom_point density try. #
    geom_smooth(
      method = "gam",
      formula = y ~ s(x, bs = "cs", k = 50),  # 'cs' gives a smooth cubic spline
      color = "#D55E00",
      linewidth = 1,
      se = TRUE,                              # keep the ribbon!
      alpha = 0.15                            # make the ribbon subtle
    ) +
    scale_x_log10(labels = comma_format()) +
    coord_cartesian(ylim = c(0, 45)) +
    labs(
      x = "Read Length (bp, log scale)",
      y = "Average Q-Score",
      title = "Read Quality vs. Read Length"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# make_quality_vs_length_plot <- function(readLengths, meanprQscore) {
#   
#   library(ggplot2)
#   library(hexbin)
#   library(patchwork)
#   
#   df <- data.frame(
#     read_length = as.numeric(readLengths),
#     avg_q_score = as.numeric(meanprQscore)
#   )
#   df <- df[df$read_length > 0 & df$avg_q_score >= 0, ]
#   
#   # Decide if log scale needed
#   use_log <- max(df$read_length) > 6000
#   y_scale <- if (use_log) scale_y_log10(labels = scales::comma) else scale_y_continuous(labels = scales::comma)
#   
#   # ---- Main hexbin ----
#   p_main <- ggplot(df, aes(x = avg_q_score, y = read_length)) +
#     stat_binhex(bins = 40) +
#     scale_fill_viridis_c(name = "Number of Reads", option = "magma") +
#     y_scale +
#     labs(x = "Phred Quality", y = "Read Length (bp)") +
#     theme_minimal(base_size = 15) +
#     theme(
#       panel.grid = element_blank(),
#       legend.position = "right"
#     )
#   
#   # ---- Top histogram (Quality) ----
#   p_top <- ggplot(df, aes(x = avg_q_score)) +
#     geom_histogram(bins = 40, fill = "grey50", color = "grey30") +
#     theme_minimal() +
#     theme(
#       axis.title.x = element_blank(),
#       axis.text.x  = element_blank(),
#       axis.ticks.x = element_blank(),
#       panel.grid = element_blank()
#     )
#   
#   # ---- Left histogram (Read length) ----
#   p_left <- ggplot(df, aes(x = read_length)) +
#     geom_histogram(bins = 40, fill = "grey50", color = "grey30") +
#     coord_flip() +
#     theme_minimal() +
#     theme(
#       axis.title.y = element_blank(),
#       axis.text.y  = element_blank(),
#       axis.ticks.y = element_blank(),
#       panel.grid = element_blank()
#     )
#   
#   # ---- Empty filler ----
#   p_empty <- ggplot() + theme_void()
#   
#   # ---- Assemble rectangular layout ----
#   final_plot <- (p_empty | p_top) /
#     (p_left  | p_main) +
#     plot_annotation(title = "Read Length vs Mean Quality") &
#     theme(plot.title = element_text(size = 18, face = "bold"))
#   
#   return(final_plot)
# }





# 6. Per-Position Quality Plot 
make_per_position_plot <- function(perPosQuality, window = 2000) {
  
  if (is.null(perPosQuality) || nrow(perPosQuality) == 0)
    return(NULL)
  
  perPosQuality$bin <- floor(perPosQuality$position / window)
  
  df <- perPosQuality |>
    dplyr::group_by(bin) |>
    dplyr::summarise(
      q10    = quantile(mean, 0.10, na.rm = TRUE),
      q25    = quantile(mean, 0.25, na.rm = TRUE),
      median = median(mean, na.rm = TRUE),
      q75    = quantile(mean, 0.75, na.rm = TRUE),
      q90    = quantile(mean, 0.90, na.rm = TRUE)
    )
  
  df$pos <- df$bin * window
  
  # Smooth curves
  smooth_loess <- function(x, y) {
    fit <- loess(y ~ x, span = 0.3)
    predict(fit)
  }
  
  df$median_s <- smooth_loess(df$pos, df$median)
  df$q25_s    <- smooth_loess(df$pos, df$q25)
  df$q75_s    <- smooth_loess(df$pos, df$q75)
  df$q10_s    <- smooth_loess(df$pos, df$q10)
  df$q90_s    <- smooth_loess(df$pos, df$q90)
  
  ggplot(df, aes(x = pos)) +
    
    # Outer 10–90% band
    geom_ribbon(
      aes(
        ymin = q10_s,
        ymax = q90_s
        # , fill = "10–90% Range"    # <-- LEGEND CODE COMMENTED
      ),
      fill = "grey70", alpha = 0.4, color = NA
    ) +
    
    # Inner 25–75% band
    geom_ribbon(
      aes(
        ymin = q25_s,
        ymax = q75_s
        # , fill = "25–75% IQR"       # <-- LEGEND CODE COMMENTED
      ),
      fill = "#56B4E9", alpha = 0.6, color = NA
    ) +
    
    # Median line
    geom_line(
      aes(
        y = median_s
        # , color = "Median (Smoothed)"  # <-- LEGEND CODE COMMENTED
      ),
      color = "#0072B2", linewidth = 1.4
    ) +
    
    labs(
      title = "Per-Read Quality Across Sequencing Read Length",
      subtitle = "Binned at 2000 bp",
      x = "Read Position (bp)",
      y = "Phred Q-Score"
    ) +
    
    theme_minimal(base_size = 15) +
    theme(
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.y = element_blank(),
      axis.line          = element_line(color = "black", linewidth = 0.6),
      axis.ticks         = element_line(color = "black", linewidth = 0.6),
      plot.title         = element_text(face = "bold", hjust = 0),
      plot.subtitle      = element_text(hjust = 0),
      legend.position    = "none"       # <-- Keep legend hidden
    ) +
    
    scale_x_continuous(
      limits = c(0, max(df$pos) + 10000),
      labels = function(x) paste0(x / 1000, "k"),
      breaks = scales::extended_breaks(n = 15),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 10),
      expand = c(0, 0)
    )
  
  # scale_fill_manual(...)   # <-- LEGEND CODE COMMENTED OUT
  # scale_color_manual(...)  # <-- LEGEND CODE COMMENTED OUT
  
  #+ #all legend related code
  # scale_fill_manual(
  #   name = "Quality Summary",
  #   values = c(
  #     "10–90% Range" = "grey70",
  #     "25–75% IQR"   = "#56B4E9"
  #   )
  # ) +
  # scale_color_manual(
  #   name = "Quality Summary",
  #   values = c(
  #     "Median (Smoothed)" = "#0072B2"
  #   )
  # )
}








# ---- Benchmark test ----
if (TRUE) {   
  t <- system.time({
    qc_obj <- QualPlot(testObj, folder_path)
  })
  
  print(t)
}


# ==== How to use: Example Usage ====

#qc_obj <- QualPlot(qc_obj, filename = "sample1.fastq")

# Visualize
# print(qc_obj@plots[["sample1.fastq"]]$length_hist)
# print(qc_obj@plots[["sample1.fastq"]]$quality_hist)
# print(qc_obj@plots[["sample1.fastq"]]$gc_hist)
#print(qc_obj@plots[["sample1.fastq"]]$q_vs_length)
# print(qc_obj@plots[["sample1.fastq"]]$per_pos_q)
