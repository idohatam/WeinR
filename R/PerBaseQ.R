#' Chunked per-position quality statistics
#'
#' Internal helper that computes per-position quality score summaries across a
#' set of reads by processing positions in chunks, which helps keep peak memory
#' use bounded for very long reads.
#'
#' @param qs_dna A Biostrings::QualityScaledDNAStringSet (or compatible) object
#'   containing per-base quality scores.
#' @param chunk_size Integer scalar giving the number of positions per chunk.
#'
#' @return A data.table with columns position, mean, median, q25, and
#'   q75.
#'
#' @keywords internal
chunked_quality_per_position <- function(qual_list, read_lengths, chunk_size = 1000) {
  
  #read_lengths <- Biostrings::width(qs_dna)
  max_len <- max(read_lengths)
  
  # Precompute order of reads by length for faster validity lookups
  ord <- order(read_lengths)
  
  chunk_stats <- vector("list", ceiling(max_len / chunk_size))
  chunk_i <- 1
  
  for (start_pos in seq(1, max_len, by = chunk_size)) {

    end_pos <- min(start_pos + chunk_size - 1, max_len)
    chunk_len <- end_pos - start_pos + 1
    
    # Fast validity lookup using monotone read lengths
    first_valid <- which(read_lengths[ord] >= start_pos)[1]
    if (is.na(first_valid)) {
      chunk_i <- chunk_i + 1
      next
    }
    
    valid_idx <- ord[first_valid:length(ord)]
    n_valid <- length(valid_idx)
    
    # Preallocate
    cm <- matrix(NA_real_, nrow = n_valid, ncol = chunk_len)
    
    # slicing
    for (row in seq_len(n_valid)) {
      i <- valid_idx[row]
      q <- qual_list[[i]]
      end_slice <- min(end_pos, read_lengths[i])
      len <- end_slice - start_pos + 1
      if (len > 0) {
        cm[row, 1:len] <- q[start_pos:end_slice]
      }
    }
    
    # One call for all quantiles
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
  
  data.table::rbindlist(chunk_stats, use.names = TRUE, fill = TRUE)
}