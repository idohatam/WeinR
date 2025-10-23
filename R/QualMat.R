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
#'  per position          TO DO
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
#'Change to non loop logic
#'use profvis on function to figure out data use stuff
#'use microbenchmark to figure out how different approaches compare

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
  quals <- quality(stringset)
  seqNames <- names(stringset)
  
  #calculate metrics
  
  #LENGTH
  lengths <- width(dna)
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
  
  #make each sequence an element in a vector -> identify longest string
  # instead of ha
  # create another vector <- nchar(vector)
  unlist vector of strings -> find largest stsring
  now we have a vector with strings and the longset value
  vector of nchar for each of those strings <- max value - that
  another vector with NA repeats
  combine into a vector will still  
  
  
    
  #Quality scores
  #per file
  qual_list <- lapply(as.character(quals), function(q) utf8ToInt(q) - 33)
  avgQscore <- mean(unlist(qual_list))
  
  #per read
  perReadQscore <- lapply(qual_list, function(q) mean(q))
                
  #try to figure out how to get it to work without NAs 
  #maybe approximate values 
  #maybe just min and max
  #why is it double
  #
  #per position 
  # pad shorter seqs with NAs <- this is really slow
  qual_mat <- do.call(rbind, lapply(qual_list,function(q){
    c(q,rep(NA, max(lengths)-length(q)))
  }))
  
  # Calculate quartiles
  meanQpb <- apply(qual_mat, 2, mean, na.rm = TRUE)
  q1Qpb <-  apply(qual_mat, 2, quantile, probs = 0.25, na.rm = TRUE)
  medQpb <-  apply(qual_mat, 2, quantile, probs = 0.50, na.rm = TRUE)
  q3Qpb <-  apply(qual_mat, 2, quantile, probs = 0.75, na.rm = TRUE)

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
  qc_obj@metrics[[file_name]] <- list(
    readLengths = lengths,
    meanprQscore = perReadQscore,
    meanppQscore = meanQpb,
    firstppQQscore = q1Qpb,
    thirdppQQscore = q3Qpb,
    medppscore = medQpb,
    Ncount = as.list(setNames(N, paste0("read", seq_along(N)))),
    prGCcontent = as.list(setNames(perReadGC, paste0("read", seq_along(perReadGC))))
  )
  
  # Add summary metrics to dataframe
  
  summary_df <- data.frame(
    file = file_name,
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
