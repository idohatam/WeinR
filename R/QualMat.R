#' calculates and fills metrics of LongReadQc object 
#'metrics: 
#'length 
#'  distribution
#'  read lengths
#'yield
#'GC content
#'  per read
#'  per position          TO DO
#'  per file    
#'N content
#'  per position          TO DO? do we even care 
#'N50
#'N90
#'Q score
#'  per read average      
#'  per position average <- takes a lot of memory, think of different approach
#'  per file average
#'  
#'  ADD ERROR HANDLING
#'  USE SYSTIME TO CHECK HOW LONG IT TAKES 
#'  FIGURE OUT WHAT IS TAKING UP SO MUCH MEMORY IN THE OBJECT
#'    save only N count if >0
#'
#'
#'use profvis on function to figure out data use stuff
#'use microbenchmark to figure out how different approaches compare
#'data table instead of dataframe?
#'try working with asciis instead of numbres?
#' try preallocating memory for memory
#' way to work not with metrics but only 
#' at x position -> take all reads with length greater than x -> calculate quartiles -> next pos
#'add adapter content as metric
#' Update workflow
#' -> try to use workflow function
#' -> generates what needs to be generated
#' -> for next week having working workflow
#'run profvis on R terminal 


QualMat <- function(qc_obj, stringset, filename){
  
  # !!! since we aren't looping, this should be done outside of the function
  # Initialize summary metrics dataframe
  #qc_obj@summary_metrics <- data.frame(
  #  file = character(),
  #  yield = numeric(),
  #  N50 = numeric(),
  #  N90 = numeric(),
  #  avgQscore = numeric(),
  #  stringsAsFactors = FALSE
  #)
  
  dna <- Biostrings::DNAStringSet(stringset)
  quals <- Biostrings::quality(stringset)
  seqNames <- names(stringset)
  
  #calculate metrics
  
  #LENGTH
  lengths <- Biostrings::width(dna)
  yield <- sum(lengths)
 
  #calculate base frequency
  af <- Biostrings::alphabetFrequency(dna)
  
  #N content per read              
  N <- af[,'N']
  
  #GC CONTENT 
  #per read
  perReadGC <- (af[,'G'] + af[,'C'])/lengths
  
  #per file
  perFileGC <- sum(af[,'G']+ af[,'C'])/yield
  
  #per position gc content
  # dataframe, 4 columns for each base: each has proportion of given base
  # row for each read
  
  
  
  #Quality scores
  #per file
  qual_list <- lapply(as.character(quals), function(q) utf8ToInt(q) - 33)
  avgQscore <- mean(unlist(qual_list))
  
  #per read
  perReadQscore <- lapply(qual_list, function(q) mean(q))
  
  # per position attempt 2 <- chunked method 
  library(matrixStats)
  
  # change chunking logic to chunk based on length.
  
  # Chunked per-position quality summary for long reads
  chunked_quality_per_position <- function(qs_dna, chunk_size = 1000) {
    qual_list <- as(qs_dna, "list")
    read_lengths <- Biostrings::width(qs_dna)
    max_len <- max(read_lengths)
    
    # Prepare list to store results per chunk
    chunk_stats <- list()
    
    # Process in chunks
    for (start_pos in seq(1, max_len, by = chunk_size)) {
      end_pos <- min(start_pos + chunk_size - 1, max_len)
      
      # Filter reads that are long enough for this chunk
      valid_idx <- which(read_lengths >= start_pos)
      if (length(valid_idx) == 0) next
      
      # Extract the relevant chunk for each valid read
      chunk_mat <- do.call(rbind, lapply(valid_idx, function(i) {
        q <- as.numeric(qual_list[[i]])
        q_chunk <- q[start_pos:min(end_pos, length(q))]
        # Pad if last chunk extends beyond read length
        if (length(q_chunk) < (end_pos - start_pos + 1)) {
          q_chunk <- c(q_chunk, rep(NA, (end_pos - start_pos + 1) - length(q_chunk)))
        }
        q_chunk
      }))
      
      # Compute per-position statistics
      mean_q <- colMeans2(chunk_mat, na.rm = TRUE)
      med_q  <- colMedians(chunk_mat, na.rm = TRUE)
      q1_q   <- colQuantiles(chunk_mat, probs = 0.25, na.rm = TRUE)
      q3_q   <- colQuantiles(chunk_mat, probs = 0.75, na.rm = TRUE)
      
      # Build results table for this chunk
      chunk_df <- data.frame(
        position = seq(start_pos, end_pos),
        mean = mean_q,
        median = med_q,
        q25 = q1_q,
        q75 = q3_q
      )
      
      chunk_stats[[length(chunk_stats) + 1]] <- chunk_df
    }
    
    # Combine all chunks
    do.call(rbind, chunk_stats)
  }
  
  # Example usage:
  q_stats <- chunked_quality_per_position(quals, chunk_size = 5000)
  # head(q_stats)
  
  #calculate summary metrics
  #N50
  decreasinglengths <- sort(lengths, decreasing = TRUE)
  tmp <- cumsum(decreasinglengths)
  N50 <- decreasinglengths[which(tmp >= yield/2)[1]]
  
  #N90
  decreasinglengths <- sort(lengths, decreasing = TRUE)
  tmp <- cumsum(decreasinglengths)
  N90 <- decreasinglengths[which(tmp >= yield*0.9)[1]]
  
  #fill object
  
  # add metrics 
  qc_obj@metrics[[filename]] <- list(
    readLengths = lengths,
    meanprQscore = perReadQscore,
    #meanppQscore = meanQpb,
    #firstppQQscore = q1Qpb,
    #thirdppQQscore = q3Qpb,
    #medppscore = medQpb,
    perPosQuality = q_stats,
    Ncount = as.list(setNames(N, paste0("read", seq_along(N)))),
    prGCcontent = as.list(setNames(perReadGC, paste0("read", seq_along(perReadGC))))
  )
  
  # Add summary metrics to dataframe
  
  summary_df <- data.frame(
    file = filename,
    yield = yield,
    N50 = N50,
    N90 = N90,
    avgQscore = avgQscore,
    stringsAsFactors = FALSE
  )
  qc_obj@summary_metrics <- rbind(qc_obj@summary_metrics, summary_df)
  
  #return the object
  base::return(qc_obj)
}
