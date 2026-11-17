library(ggplot2)
library(Rsamtools)

#' Process per-read modification density
#'
#' Calculates recorded modification calls per read as a percentage of read length.
#'
#' @param bam_file Path to the BAM file.
#'
#' @return Numeric vector of per-read modification densities.
#' @keywords internal
process_modification_density <- function(bam_file) {
  
  bam_data <- Rsamtools::scanBam(
    bam_file,
    param = Rsamtools::ScanBamParam(
      what = "seq",
      tag = "MM"
    )
  )
  
  read_lengths <- nchar(as.character(bam_data[[1]]$seq))
  mm_tags <- bam_data[[1]]$tag$MM
  
  modification_counts <- numeric(length(mm_tags))
  
  for (i in seq_along(mm_tags)) {
    
    if (is.na(mm_tags[i]) || mm_tags[i] == "") {
      next
    }
    
    groups <- strsplit(mm_tags[i], ";", fixed = TRUE)[[1]]
    
    for (group in groups) {
      
      if (group == "") {
        next
      }
      
      parts <- strsplit(group, ",", fixed = TRUE)[[1]]
      
      modification_counts[i] <- modification_counts[i] +
        length(parts) - 1
    }
  }
  
  read_density <- modification_counts / read_lengths * 100
  read_density[!is.finite(read_density)] <- NA_real_
  
  read_density
}


#' Overall base modification density histogram
#'
#' Produces a histogram of recorded modification calls per read, expressed as
#' a percentage of read length.
#'
#' @param read_density Numeric vector of per-read modification densities.
#' @param file_name Character title for the plot.
#'
#' @return A ggplot object, or NULL if no finite values are available.
#' @keywords internal
make_mod_density_plot <- function(read_density, file_name) {
  
  read_density <- as.numeric(base::unlist(read_density, use.names = FALSE))
  read_density <- read_density[base::is.finite(read_density)]
  
  if (base::length(read_density) == 0L) {
    return(NULL)
  }
  
  ggplot2::ggplot(
    data = data.frame(modification_density = read_density),
    ggplot2::aes(x = modification_density)
  ) +
    ggplot2::geom_histogram(
      bins = 50,
      fill = "#0072B2",
      color = "white",
      alpha = 0.9
    ) +
    ggplot2::geom_text(
      stat = "bin",
      ggplot2::aes(
        label = ifelse(
          ggplot2::after_stat(count) > 0,
          ggplot2::after_stat(count),
          ""
        )
      ),
      bins = 50,
      vjust = -0.3,
      size = 1.8
    ) +
    ggplot2::scale_x_continuous(
      name = "Recorded Modification Calls per Read (%)",
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 10)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 10)
    ) +
    ggplot2::labs(
      title = file_name,
      y = "Read Count"
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.6),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      legend.position = "none",
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(12, 15, 15, 15)
    )
}


# bam_file <- "/Users/aaronsantos/Downloads/subsampled_1percent.bam"
# 
# read_density <- process_modification_density(bam_file)
# 
# print(
#   make_mod_density_plot(
#     read_density,
#     "PBM60157_1pct Base Modification Density per Read"
#   )
# )