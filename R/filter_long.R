

filter_long <- function(chunk_qsds, single_qc, filter_summary, start_index, end_index, min_avg_q, min_length, max_n) {
  
  filename <- single_qc@files[[1]]
  key <- basename(filename)
  
  if (length(single_qc@files) != 1L) {
    stop("FilterLong() expects a QC object with exactly one file in qc_obj@files.", call. = FALSE)
  }
  file_metrics <- single_qc@metrics[[filename]]
  if (is.null(file_metrics)) {
    stop("Metrics not found for file: ", filename, call. = FALSE)
  }
  
  lengths  <- (file_metrics$readLengths)[start_index:end_index]
  q_scores <- (unlist(file_metrics$meanprQscore))[start_index:end_index]
  n_counts <- (unlist(file_metrics$Ncount))[start_index:end_index]
  
  filter <- (q_scores >= min_avg_q) &
    (lengths  >= min_length) &
    (n_counts <= max_n)
  
  result <- chunk_qsds[filter]
  filter_summary$reads_before <- filter_summary$reads_before + length(chunk_qsds)
  filter_summary$reads_after <- filter_summary$reads_after + length(result)
  
  processed_chunk <- result
}