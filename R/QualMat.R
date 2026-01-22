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
  
  seqNames <- names(stringset)
  
  #calculate metrics
  
  #LENGTH
  lengths <- Biostrings::width(Biostrings::DNAStringSet(stringset))
  yield <- sum(lengths)
 
  #calculate base frequency
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
  qual_list <- lapply(as.character(Biostrings::quality(stringset)), function(q) utf8ToInt(q) - 33)
  avgQscore <- mean(unlist(qual_list))
  
  #per read
  perReadQscore <- lapply(qual_list, function(q) mean(q))
  
  #per position
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
  
  # Add summary metrics 
  summary_df <- data.frame(
    file = filename,
    yield = yield,
    N50 = N50,
    N90 = N90,
    avgQscore = avgQscore,
    `N count` = sum(N),
    stringsAsFactors = FALSE
  )
  qc_obj@summary_metrics <- rbind(qc_obj@summary_metrics, summary_df)
  
  #return the object
  base::return(qc_obj)
}

library(parallel)
                                                 # what does the chunk mean????
chunked_quality_per_position <- function(qs_dna, chunk_size = 1000) {
  
  print('convert to list')
  start <- Sys.time()
  # Convert once
  #qual_list <- lapply(as(qs_dna, "list"), as.numeric)
 
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, "qs_list")  # export your object
  clusterEvalQ(cl, library(Biostrings))  # make sure Biostrings is loaded
  
  qual_list <- parLapply(cl, seq_along(qs_dna), function(i) {
    as.numeric(qs_dna[[i]])
  })
  
  stopCluster(cl)
  
  end <- Sys.time()
  print(end - start)
  
  print('calculate max length')
  start <- Sys.time()
  read_lengths <- Biostrings::width(qs_dna)
  max_len <- max(read_lengths)
  end <- Sys.time()
  print(end - start)
  
  print('order reads ')
  start <- Sys.time()
  # Precompute order of reads by length for faster validity lookups
  ord <- order(read_lengths)
  end <- Sys.time()
  print(end - start)
  
  chunk_stats <- vector("list", ceiling(max_len / chunk_size))
  chunk_i <- 1
  
  loopstart <- Sys.time()
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
    
    #slicing
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
    quants <- matrixStats::colQuantiles(cm, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    q25 <- quants[, 1]
    med <- quants[, 2]
    q75 <- quants[, 3]
    
    # mean still needs its own pass
    meanv <- matrixStats::colMeans2(cm, na.rm = TRUE)
    
    chunk_stats[[chunk_i]] <- data.table::data.table(
      position = start_pos:end_pos,
      mean     = meanv,
      median   = med,
      q25      = q25,
      q75      = q75
    )
    
    chunk_i <- chunk_i + 1
  }
  loopend <- Sys.time()
  print('time looping ')
  print(loopend - loopstart)
  
  print('binding')
  start <- Sys.time()
  data.table::rbindlist(chunk_stats)
  end <- Sys.time()
  print(end - start)
}
