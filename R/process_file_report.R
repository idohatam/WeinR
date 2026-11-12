test <- "hello"
print(test)

process_file_report <- function(filename, qc_obj) {
  
  sourceCpp("./C/weinr.cpp")
  
  CHUNK_SIZE <- 10000L
  MAX_WIDTH <- 2000000L
  BIN_WIDTH <- 5
  BINS <- 20
  
  isBAM <- FALSE
  per_read_lengths <- c()
  per_read_N_content <- c()
  per_read_GC_content <- c()
  per_read_avg_error_p <- c()
  per_read_avg_q_score <- c()
  chunk_qpp_subtotals <- vector(mode = "list")
  qpp_count <- integer(MAX_WIDTH)
  qpp_sum <- numeric(MAX_WIDTH)
  qpp_frequency_table <- matrix(0, nrow = BINS, ncol = MAX_WIDTH)
  
  if (grepl("\\.bam$", filename, ignore.case = TRUE)) {
    bam_file <- filename
    bam_stream <- open(BamFile(bam_file, yieldSize = CHUNK_SIZE))
    isBAM <- TRUE
  } else {
    fastq_file <- filename
    fastq_stream <- FastqStreamer(fastq_file, n = CHUNK_SIZE)
    isBAM <- FALSE
  }
  
  start_job <- Sys.time()
  start_time <- start_job
  
  repeat {
    
    if(isBAM) {
      chunk <- scanBam(bam_stream)
      if(length(chunk[[1]][[1]]) == 0) {
        break
      }
      chunk_qsds <- QualityScaledDNAStringSet(chunk[[1]]$seq, chunk[[1]]$qual)
    } else {
      chunk <- yield(fastq_stream)
      if(length(chunk) == 0) {
        break
      }
      chunk_qsds <- as(chunk, "QualityScaledDNAStringSet")
    }
    
    # AGGREGATE DNA METRICS
    
    chunk_dna_ss <- Biostrings::DNAStringSet(chunk_qsds)
    
    chunk_per_read_lengths <- Biostrings::width(chunk_dna_ss)
    per_read_lengths = c(per_read_lengths, chunk_per_read_lengths)
    
    chunk_per_read_base_frequency <- Biostrings::alphabetFrequency(chunk_dna_ss)
    
    chunk_per_read_N_content <- chunk_per_read_base_frequency[,'N']
    per_read_N_content <- c(per_read_N_content, chunk_per_read_N_content)
    
    chunk_per_read_G_content <- chunk_per_read_base_frequency[,'G']
    chunk_per_read_C_content <- chunk_per_read_base_frequency[,'C']
    chunk_per_read_GC_content <- (chunk_per_read_G_content + chunk_per_read_C_content)/chunk_per_read_lengths #reported as percentage
    per_read_GC_content <- c(per_read_GC_content, chunk_per_read_GC_content)
    
    # AGGREGATE QUALITY METRICS
    
    chunk_quality_ss <- as.character(Biostrings::quality(chunk_qsds))
    chunk_error_p <- sanger_to_probability(chunk_quality_ss)
    
    chunk_per_read_avg_q_score <- calculate_per_read_avg_q_score(chunk_error_p)
    per_read_avg_q_score <- c(per_read_avg_q_score, chunk_per_read_avg_q_score)
    
    chunk_qpp_subtotals <- tabulate_qpp_mean(chunk_error_p, MAX_WIDTH)
    qpp_count = qpp_count + chunk_qpp_subtotals[[1]]
    qpp_sum = qpp_sum + chunk_qpp_subtotals[[2]]
    
    qpp_frequency_table <- tabulate_qpp_median(qpp_frequency_table, chunk_quality_ss, BIN_WIDTH)
    
  }
  
  # CALCULATE FILE METRICS
  
  yield <- sum(per_read_lengths)
  read_count <- length(per_read_lengths)
  N_count <- sum(per_read_N_content)
  decreasinglengths <- sort(per_read_lengths, decreasing = TRUE)
  tmp <- cumsum(as.numeric(decreasinglengths))
  N50 <- decreasinglengths[which(tmp >= yield/2)[1]]
  N90 <- decreasinglengths[which(tmp >= yield*0.9)[1]]
  
  avg_Q_score <- -10 * log10(sum(qpp_sum) / sum(qpp_count))
  qpp_mean <- calculate_qpp_mean(qpp_count, qpp_sum, MAX_WIDTH)
  qpp_median <- calculate_qpp_median(qpp_frequency_table, MAX_WIDTH, BINS, BIN_WIDTH)
  
  
  if(isBAM) {
    close(bam_stream)
  } else {
    close(fastq_stream)
  }
  
  end_job <- Sys.time()
  print(end_job - start_job)
  
  # STORE METRICS IN qc_obj
  
  MAX_LENGTH = max(per_read_lengths)
  mean = qpp_mean[1:MAX_LENGTH]
  median = qpp_median[1:MAX_LENGTH]
  
  q_stats <- data.table(
    position = 1:MAX_LENGTH,
    mean = mean,
    median = median,
    q25 = double(MAX_LENGTH),
    q75 = double(MAX_LENGTH)
  )
  
  qc_obj@metrics[[filename]] <- list(
    readLengths = per_read_lengths,
    meanprQscore = as.list(per_read_avg_q_score),
    perPosQuality = q_stats,
    Ncount = as.list(per_read_N_content),
    prGCcontent = as.list(per_read_GC_content)
  )
  
  summary_df <- data.frame(
    File = filename,
    `Read Count` = read_count,
    `Yield (bp)` = yield,
    N50 = N50,
    N90 = N90,
    `Average Per Read Q Score` = round(avg_Q_score, digits = 2),
    `N Count` = N_count,
    stringsAsFactors = FALSE
  )
  qc_obj@summary_metrics <- rbind(qc_obj@summary_metrics, summary_df)
  
  qc_obj
}