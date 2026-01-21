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
  
  #dna <- Biostrings::DNAStringSet(stringset)
  #quals <- Biostrings::quality(stringset)
  seqNames <- names(stringset)
  
  #calculate metrics
  
  #LENGTH
  #lengths <- Biostrings::width(dna)
  #yield <- sum(lengths)
  lengths <- Biostrings::width(Biostrings::DNAStringSet(stringset))
  yield <- sum(lengths)
 
  #calculate base frequency
  #af <- Biostrings::alphabetFrequency(dna)
  af <- Biostrings::alphabetFrequency(Biostrings::DNAStringSet(stringset))
  
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
  #qual_list <- lapply(as.character(quals), function(q) utf8ToInt(q) - 33)
  qual_list <- lapply(as.character(Biostrings::quality(stringset)), function(q) utf8ToInt(q) - 33)
  avgQscore <- mean(unlist(qual_list))
  
  #per read
  perReadQscore <- lapply(qual_list, function(q) mean(q))
  
  #calculate per position messages
  q_stats <- chunked_quality_per_position(Biostrings::quality(stringset), chunk_size = 5000)
  
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
    perPosQuality = q_stats,
    Ncount = as.list(N),
    prGCcontent = as.list(perReadGC)
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
      mean_q <- matrixStats::colMeans2(chunk_mat, na.rm = TRUE)
      med_q  <- matrixStats::colMedians(chunk_mat, na.rm = TRUE)
      q1_q   <- matrixStats::colQuantiles(chunk_mat, probs = 0.25, na.rm = TRUE)
      q3_q   <- matrixStats::colQuantiles(chunk_mat, probs = 0.75, na.rm = TRUE)
      
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


BenchmarkQualMat <- function(qc_obj, stringset, filename){
  functionstart <- Sys.time()
  #dna <- Biostrings::DNAStringSet(stringset)
  #quals <- Biostrings::quality(stringset)
  seqNames <- names(stringset)
  
  #calculate metrics
  
  #LENGTH
  #lengths <- Biostrings::width(dna)
  #yield <- sum(lengths)
  lt <- Sys.time()
  lengths <- Biostrings::width(Biostrings::DNAStringSet(stringset))
  yield <- sum(lengths)
  lte <- Sys.time()
  print('length and yield calculation:')
  print(lte - lt)
  
  #calculate base frequency
  #af <- Biostrings::alphabetFrequency(dna)
  aft <- Sys.time()
  af <- Biostrings::alphabetFrequency(Biostrings::DNAStringSet(stringset))
  afte <- Sys.time()
  print('alphabetfrequency calculations:')
  print(afte - aft)
  
  
  #N content per read  
  ns <- Sys.time()
  N <- af[,'N']
  
  #GC CONTENT 
  #per read
  perReadGC <- (af[,'G'] + af[,'C'])/lengths
  
  #per file
  perFileGC <- sum(af[,'G']+ af[,'C'])/yield
  afe <- Sys.time()
  print('N and GC calculations:')
  print(afe - ns)
  #per position gc content
  # dataframe, 4 columns for each base: each has proportion of given base
  # row for each read
  
  #Quality scores
  #per file
  #qual_list <- lapply(as.character(quals), function(q) utf8ToInt(q) - 33)
  qt <- Sys.time()
  qual_list <- lapply(as.character(Biostrings::quality(stringset)), function(q) utf8ToInt(q) - 33)
  qle <- Sys.time()
  print('making qual list')
  print(qle - qt)
  
  ms <- Sys.time()
  avgQscore <- mean(unlist(qual_list))
  me <- Sys.time()
  print('per file quality scores')
  print(me - ms)
  
  #per read
  qre <- Sys.time()
  perReadQscore <- lapply(qual_list, function(q) mean(q))
  qe <- Sys.time()
  print('per read quality scores')
  print(qe - qre)
  
  #calculate per position messages
  qs <- Sys.time()
  q_stats <- chunked_quality_per_position(Biostrings::quality(stringset), chunk_size = 5000)
  qse <- Sys.time()
  print('per position q score calculation:')
  print(qse - qs)
  
  #calculate summary metrics
  #N50
  ns <- Sys.time()
  decreasinglengths <- sort(lengths, decreasing = TRUE)
  tmp <- cumsum(decreasinglengths)
  N50 <- decreasinglengths[which(tmp >= yield/2)[1]]
  
  #N90
  decreasinglengths <- sort(lengths, decreasing = TRUE)
  tmp <- cumsum(decreasinglengths)
  N90 <- decreasinglengths[which(tmp >= yield*0.9)[1]]
  ne <- Sys.time()
  print('N50/90 calculations')
  print(ne - ns)
  
  #fill object
  # add metrics 
  fs <- Sys.time()
  qc_obj@metrics[[filename]] <- list(
    readLengths = lengths,
    meanprQscore = perReadQscore,
    perPosQuality = q_stats,
    Ncount = as.list(N),
    prGCcontent = as.list(perReadGC)
  )
  fe <- Sys.time()
  print('fill metrics')
  print(fe-fs)
  
  # Add summary metrics to dataframe
  ss <- Sys.time()
  summary_df <- data.frame(
    file = filename,
    yield = yield,
    N50 = N50,
    N90 = N90,
    avgQscore = avgQscore,
    stringsAsFactors = FALSE
  )
  qc_obj@summary_metrics <- rbind(qc_obj@summary_metrics, summary_df)
  se <- Sys.time()
  print('fill summary metrics:')
  print(se - ss)
  
  #return the object
  fe <- Sys.time()
  print('overall runtime')
  print(fe - functionstart) 
  base::return(qc_obj)
}

chunked_quality_per_position <- function(qs_dna, chunk_size = 1000) {
  
  # Convert once
  qual_list <- lapply(as(qs_dna, "list"), as.numeric)
  read_lengths <- Biostrings::width(qs_dna)
  max_len <- max(read_lengths)
  
  # Precompute order of reads by length for faster validity lookups
  ord <- order(read_lengths)
  
  chunk_stats <- vector("list", ceiling(max_len / chunk_size))
  chunk_i <- 1
  
  for (start_pos in seq(1, max_len, by = chunk_size)) {
    
    end_pos <- min(start_pos + chunk_size - 1, max_len)
    chunk_len <- end_pos - start_pos + 1
    
    # Fast validity lookup using monotone read lengths
    # Reads are sorted so we only find the cutoff once
    first_valid <- which(read_lengths[ord] >= start_pos)[1]
    if (is.na(first_valid)) {
      chunk_i <- chunk_i + 1
      next
    }
    
    valid_idx <- ord[first_valid:length(ord)]
    n_valid <- length(valid_idx)
    
    # Preallocate
    cm <- matrix(NA_real_, nrow = n_valid, ncol = chunk_len)
    
    # 🔥 Fast slicing
    for (row in seq_len(n_valid)) {
      i <- valid_idx[row]
      q <- qual_list[[i]]
      end_slice <- min(end_pos, read_lengths[i])
      len <- end_slice - start_pos + 1
      if (len > 0) {
        cm[row, 1:len] <- q[start_pos:end_slice]
      }
    }
    
    #One call for all quantiles
    quants <- colQuantiles(cm, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    q25 <- quants[, 1]
    med <- quants[, 2]
    q75 <- quants[, 3]
    
    # mean still needs its own pass
    meanv <- colMeans2(cm, na.rm = TRUE)
    
    chunk_stats[[chunk_i]] <- data.table(
      position = start_pos:end_pos,
      mean     = meanv,
      median   = med,
      q25      = q25,
      q75      = q75
    )
    
    chunk_i <- chunk_i + 1
  }
  
  rbindlist(chunk_stats)
}

