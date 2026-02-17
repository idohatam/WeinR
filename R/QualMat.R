#' Compute and store per-file quality-control metrics
#'
#' Internal helper that computes read-level and file-level summary metrics from a
#' `DNAStringSet`/`QualityScaledDNAStringSet`-like input and stores them into a
#' QC object.
#'
#' @param qc_obj An S4 QC object with at least `@metrics` and `@summary_metrics`
#'   slots. This function appends metrics for `filename` into these slots.
#' @param stringset A `Biostrings` string set containing the sequences (and
#'   optionally qualities, depending on object type).
#' @param filename Character scalar naming the file being summarized (used as a
#'   key in `qc_obj@metrics` and in the summary table).
#'
#' @return The updated `qc_obj`.
#'
#' @keywords internal
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
  
  rm(af)
  
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
  q_stats <- chunked_quality_per_position(qual_list,lengths, chunk_size = 1000)
  
  rm(qual_list)
  
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
