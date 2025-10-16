#' calculates and fills metrics of LongReadQc object 
#'metrics: 
#'length 
#'yield
#'GC content
#'  per read
#'  per position  TO DO
#'  per file      TO DO
#'N content
#'N50
#'N90
#'Q score
#'  per read average      TODO
#'  per position average 
#'  per file average

# Load required package
library(Biostrings)

QualMat <- function(qc_obj, listofstringsets){
  
  # Initialize summary metrics dataframe
  summary_df <- data.frame(
    file = character(),
    yield = numeric(),
    N50 = numeric(),
    N90 = numeric(),
    avgQscore = numeric(),
    stringsAsFactors = FALSE
  )
  
  #iterate through the stringset
  for(i in seq_along(listofstringsets)){
    
    stringset <- listofstringsets[[i]]
    
    #get file names 
    file_name <- qc_obj@files[i]
    
    dna <- DNAStringSet(stringset)
    quals <- quality(stringset)
    seqNames <- names(stringset)
    
    #calculate metrics
    
    #LENGTH
    lengths <- width(dna)
    yield <- sum(lengths)
   
     
    #calculate base frequency
    af <- alphabetFrequency(dna)
    
    #N content per read              
    N <- af[,'N']
    
    #GC CONTENT 
    #per read
    perReadGC <- (af[,'G'] + af[,'C'])/lengths
    
    #per file
    
    #per position gc content
    # dataframe, 4 columns for each base: each has proportion of given base
    # row for each read
    
    #Quality scores
    #per file
    qual_list <- lapply(as.character(quals), function(q) utf8ToInt(q) - 33)
    avgQscore <- mean(unlist(qual_list))
    
    #per read
                  
    #per position 
    # pad shorter seqs with NAs
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
      ppQscore = meanQpb,
      Ncount = as.list(setNames(N, paste0("read", seq_along(N)))),
      prGCcontent = as.list(setNames(perReadGC, paste0("read", seq_along(perReadGC))))
    )
    
    # Add summary metrics to dataframe
    summary_df <- rbind(summary_df, data.frame(
      file = file_name,
      yield = yield,
      N50 = N50,
      N90 = N90,
      avgQscore = avgQscore,
      stringsAsFactors = FALSE
    ))
  }
  #return the object
  base::return(qc_obj)
}


