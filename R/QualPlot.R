# File-by-File Mode: QualPlot() Definition
# Compatible with LongReadQC structure

# Required libraries
library(ggplot2)
library(scales)
library(mgcv) # for geom_smooth(method = "gam")
library(ggrastr)

# Main function: process one file at a time
QualPlot <- function(qc_obj, filename) {
  if (!filename %in% names(qc_obj@metrics)) {
    stop(paste("File", filename, "not found in qc_obj@metrics"))
  }
  
  # Extract metrics list for this file
  metrics <- qc_obj@metrics[[filename]]
  
  # Generate and store all plots (per file)
  qc_obj@plots[[filename]] <- list(
    length_hist   = make_length_plot(metrics$readLengths),
    quality_hist  = make_quality_plot(metrics$meanprQscore),
    gc_hist       = make_gc_plot(metrics$prGCcontent),
    q_vs_length   = make_quality_vs_length_plot(metrics$readLengths, metrics$meanprQscore),
    per_pos_q     = make_per_position_plot(metrics$perPosQuality)
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

# 5. Q Score vs. Read Length (Scatter + Smooth)
make_quality_vs_length_plot <- function(readLengths, meanprQscore) {
  df <- data.frame(
    read_length = as.numeric(unlist(readLengths, use.names = FALSE)),
    avg_q_score = as.numeric(unlist(meanprQscore, use.names = FALSE))
  )
  
  ggplot(df, aes(x = read_length, y = avg_q_score)) +
    geom_point_rast(alpha = 0.25, color = "#0072B2", size = 0.4) +
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

# 6. Per-Position Quality Plot
make_per_position_plot <- function(perPosQuality, window = 100, save_path = NULL) {
  if (is.null(perPosQuality) || nrow(perPosQuality) == 0) return(NULL)
  
  perPosQuality <- perPosQuality[order(perPosQuality$position), ]
  perPosQuality$bin <- floor(perPosQuality$position / window)
  df <- aggregate(
    cbind(mean, q25, q75) ~ bin,
    data = perPosQuality,
    FUN = mean,
    na.rm = TRUE
  )
  names(df) <- c("position_bin", "mean", "q25", "q75")
  
  p <- ggplot(df, aes(x = position_bin)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightblue", alpha = 0.4) +
    geom_line(aes(y = mean), color = "blue", linewidth = 0.4) +
    labs(
      x = sprintf("Read Position (binned every %d bases)", window),
      y = "Quality Score",
      title = "Per-Position Quality"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  return(p)
}

# ==== How to use: Example Usage ====

#qc_obj <- QualPlot(qc_obj, filename = "sample1.fastq")

# Visualize
# print(qc_obj@plots[["sample1.fastq"]]$length_hist)
# print(qc_obj@plots[["sample1.fastq"]]$quality_hist)
# print(qc_obj@plots[["sample1.fastq"]]$gc_hist)
#print(qc_obj@plots[["sample1.fastq"]]$q_vs_length)
# print(qc_obj@plots[["sample1.fastq"]]$per_pos_q)
