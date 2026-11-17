library(ggplot2)
library(Rsamtools)

#' Process per-read modification probabilities
#'
#' Extracts ML probabilities for each modification type and calculates
#' per-read mean probabilities and file-wide mean probabilities.
#'
#' @param bam_file Path to the BAM file.
#'
#' @return A list containing modification probabilities, per-read means,
#'   and file-wide means.
#' @keywords internal
process_modification_probabilities <- function(bam_file) {
  
  bam_data <- Rsamtools::scanBam(
    bam_file,
    param = Rsamtools::ScanBamParam(
      tag = c("MM", "ML")
    )
  )
  
  mm_tags <- bam_data[[1]]$tag$MM
  ml_tags <- bam_data[[1]]$tag$ML
  
  modifications <- c(
    "m5C" = "^C\\+m",
    "m6A" = "^A\\+a",
    "inosine" = "^A\\+17596",
    "pseU" = "^T\\+17802",
    "2OmeC" = "^C\\+19228",
    "2OmeA" = "^A\\+69426",
    "2OmeG" = "^G\\+19229",
    "2OmeU" = "^T\\+19227"
  )
  
  type_values <- setNames(
    vector("list", length(modifications)),
    names(modifications)
  )
  
  read_means <- setNames(
    vector("list", length(modifications)),
    names(modifications)
  )
  
  for (i in seq_along(mm_tags)) {
    
    mm <- mm_tags[[i]]
    ml <- as.numeric(ml_tags[[i]])
    
    if (length(mm) == 0 || is.na(mm[1]) || mm[1] == "") {
      next
    }
    
    groups <- strsplit(mm[1], ";", fixed = TRUE)[[1]]
    
    ml_position <- 1
    
    read_type_values <- setNames(
      vector("list", length(modifications)),
      names(modifications)
    )
    
    for (group in groups) {
      
      if (group == "") {
        next
      }
      
      parts <- strsplit(group, ",", fixed = TRUE)[[1]]
      
      modification_code <- parts[1]
      number_of_calls <- length(parts) - 1
      
      if (number_of_calls == 0) {
        next
      }
      
      if (ml_position + number_of_calls - 1 > length(ml)) {
        break
      }
      
      probabilities <- ml[
        ml_position:(ml_position + number_of_calls - 1)
      ]
      
      ml_position <- ml_position + number_of_calls
      
      modification_type <- names(modifications)[
        vapply(
          modifications,
          function(pattern) grepl(pattern, modification_code),
          logical(1)
        )
      ]
      
      if (length(modification_type) > 0) {
        
        modification_type <- modification_type[1]
        
        read_type_values[[modification_type]] <- c(
          read_type_values[[modification_type]],
          probabilities
        )
      }
    }
    
    for (modification_type in names(modifications)) {
      
      values <- read_type_values[[modification_type]]
      
      if (length(values) > 0) {
        
        type_values[[modification_type]] <- c(
          type_values[[modification_type]],
          values
        )
        
        read_means[[modification_type]] <- c(
          read_means[[modification_type]],
          mean(values / 255)
        )
      }
    }
  }
  
  average_ML_file <- lapply(
    type_values,
    function(values) {
      if (length(values) > 0) {
        mean(values / 255)
      } else {
        NA_real_
      }
    }
  )
  
  list(
    type_values = type_values,
    read_means = read_means,
    average_ML_file = average_ML_file
  )
}


#' Modification probability histogram
#'
#' Produces a histogram of per-read mean modification probabilities.
#'
#' @param average_ML_read Numeric vector of per-read mean probabilities.
#' @param modification_type Modification type label.
#'
#' @return A ggplot object, or NULL if no valid probabilities are available.
#' @keywords internal
make_modification_probability_plot <- function(average_ML_read, modification_type) {
  
  average_ML_read <- as.numeric(
    base::unlist(average_ML_read, use.names = FALSE)
  )
  
  average_ML_read <- average_ML_read[
    base::is.finite(average_ML_read)
  ]
  
  if (base::length(average_ML_read) == 0L) {
    return(NULL)
  }
  
  ggplot2::ggplot(
    data = data.frame(average_ML = average_ML_read),
    ggplot2::aes(x = average_ML)
  ) +
    ggplot2::geom_histogram(
      bins = 50,
      fill = "#0072B2",
      color = "white",
      alpha = 0.9
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(),
      breaks = scales::pretty_breaks(n = 10),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 8)
    ) +
    ggplot2::labs(
      x = paste("Average", modification_type, "Probability per Read"),
      y = "Number of Reads",
      title = paste(
        modification_type,
        "Modification Probability Distribution"
      )
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0
      ),
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = 0.6
      ),
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = 0.6
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(12, 15, 15, 15)
    )
}


# bam_file <- "/Users/aaronsantos/Downloads/subsampled_1percent.bam"
# 
# modification_data <- process_modification_probabilities(bam_file)
# 
# for (modification_type in names(modification_data$type_values)) {
#   
#   if (length(modification_data$type_values[[modification_type]]) > 0) {
#     
#     print(modification_type)
#     print(modification_data$average_ML_file[[modification_type]])
#     
#     print(
#       make_modification_probability_plot(
#         modification_data$read_means[[modification_type]],
#         modification_type
#       )
#     )
#   }
# }