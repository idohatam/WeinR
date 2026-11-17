library(ggplot2)
library(Rsamtools)

#' Count base frequencies
#'
#' Counts the frequency of each nucleotide in the BAM file.
#'
#' @param bam_file Path to the BAM file.
#'
#' @return Data frame containing nucleotide frequencies.
#' @keywords internal
process_base_frequencies <- function(bam_file) {
  
  bam_data <- scanBam(
    bam_file,
    param = ScanBamParam(
      what = "seq"
    )
  )
  
  sequences <- as.character(bam_data[[1]]$seq)
  sequences <- toupper(paste(sequences, collapse = ""))
  
  bases <- strsplit(sequences, "")[[1]]
  base_frequencies <- table(factor(bases, levels = c("A", "C", "G", "T", "N")))
  
  data.frame(
    Base = names(base_frequencies),
    Frequency = as.numeric(base_frequencies)
  )
}


#' Plot base frequencies
#'
#' Produces a bar plot of nucleotide frequencies.
#'
#' @param base_frequencies Data frame of nucleotide frequencies.
#' @param file_name Character title for the plot.
#'
#' @return A ggplot object.
#' @keywords internal
make_base_frequency_plot <- function(base_frequencies, file_name) {
  
  ggplot(
    base_frequencies,
    aes(x = Base, y = Frequency)
  ) +
    geom_col(
      fill = "#0072B2",
      color = "white",
      alpha = 0.9
    ) +
    scale_y_continuous(
      labels = scales::comma_format(),
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 12)
    ) +
    labs(
      x = "Base",
      y = "Frequency",
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


# bam_file <- "/Users/aaronsantos/Downloads/subsampled_1percent.bam"
# 
# base_frequencies <- process_base_frequencies(bam_file)
# 
# print(
#   make_base_frequency_plot(
#     base_frequencies,
#     "Base Frequency Distribution"
#   )
# )