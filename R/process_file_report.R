# ---- BAM modification settings (used by parse_chunk and process_file_report) ----

KEEP_ML_VALUES <- TRUE    # FALSE = don't store every raw ML value (saves lots of RAM);
# average_ML_file is computed from running sums either way

modifications <- c(
  "m5C"     = "^C\\+m",
  "m6A"     = "^A\\+a",
  "inosine" = "^A\\+17596",
  "pseU"    = "^T\\+17802",
  "2OmeC"   = "^C\\+19228",
  "2OmeA"   = "^A\\+69426",
  "2OmeG"   = "^G\\+19229",
  "2OmeU"   = "^T\\+19227"
)
mod_names <- names(modifications)
n_mod     <- length(modifications)

# ---- helpers ----------------------------------------------------------------

# scanBam returns tags as a vector, or a list with NULL for reads lacking the tag
tag_first <- function(x, na) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    vapply(x, function(v) if (length(v) == 0L) na else v[[1L]], FUN.VALUE = na)
  } else x
}

# Sum `vals` by integer group id `g` into a length-n vector
sum_by_read <- function(vals, g, n) {
  out <- numeric(n)
  if (length(g) > 0L) {
    agg <- rowsum(vals, g, reorder = FALSE)
    out[as.integer(rownames(agg))] <- agg[, 1L]
  }
  out
}

# Parse MM/ML for a whole chunk at once (no per-read, per-modification loops).
# Semantics match the original: groups with no calls are skipped, and groups whose
# calls would run past the end of the read's ML array are dropped.
parse_chunk <- function(mm_tags, ml_tags, read_lengths) {
  n <- length(read_lengths)
  if (is.null(mm_tags)) mm_tags <- rep(NA_character_, n)
  if (is.null(ml_tags)) ml_tags <- vector("list", n)
  
  mm      <- tag_first(mm_tags, NA_character_)
  ml_len  <- lengths(ml_tags)
  ml_all  <- as.numeric(unlist(ml_tags, use.names = FALSE))
  ml_off  <- c(0, cumsum(ml_len))[seq_len(n)]      # 0-based start of each read in ml_all
  cs      <- c(0, cumsum(ml_all))                  # prefix sums for O(1) group sums
  
  has_mm  <- which(!is.na(mm) & nzchar(mm))
  groups  <- strsplit(mm[has_mm], ";", fixed = TRUE)
  g       <- unlist(groups, use.names = FALSE)
  g_read  <- rep(has_mm, lengths(groups))
  
  ncalls  <- nchar(g) - nchar(gsub(",", "", g, fixed = TRUE))
  keep    <- nzchar(g) & ncalls > 0
  g <- g[keep]; g_read <- g_read[keep]; ncalls <- ncalls[keep]
  
  empty <- list(
    density = lapply(seq_len(n_mod), function(i) rep(NA_real_, n)),
    means   = vector("list", n_mod),
    values  = vector("list", n_mod),
    sum     = numeric(n_mod),
    n       = numeric(n_mod)
  )
  names(empty$density) <- names(empty$means) <- names(empty$values) <- mod_names
  names(empty$sum) <- names(empty$n) <- mod_names
  if (length(g) == 0L) {
    return(empty)
  }
  
  # position of each group inside its read's ML array
  code    <- sub(",.*", "", g)
  cum     <- cumsum(ncalls)
  r       <- rle(g_read)
  first   <- cumsum(r$lengths) - r$lengths + 1L
  base    <- rep(cum[first] - ncalls[first], r$lengths)
  end_in  <- cum - base
  start_in <- end_in - ncalls
  valid   <- end_in <= ml_len[g_read]
  
  # classify each *unique* mod code once (first matching pattern wins)
  ucodes  <- unique(code)
  utype   <- rep(NA_integer_, length(ucodes))
  for (k in rev(seq_len(n_mod))) utype[grepl(modifications[[k]], ucodes)] <- k
  type_idx <- utype[match(code, ucodes)]
  
  sel <- valid & !is.na(type_idx)
  g_read <- g_read[sel]; ncalls <- ncalls[sel]; type_idx <- type_idx[sel]
  gstart <- ml_off[g_read] + start_in[sel] + 1          # 1-based into ml_all
  gsum   <- cs[gstart + ncalls] - cs[gstart]            # sum of ML probs per group
  
  out <- empty
  for (k in seq_len(n_mod)) {
    i <- which(type_idx == k)
    cnt <- sum_by_read(ncalls[i], g_read[i], n)
    sm  <- sum_by_read(gsum[i],   g_read[i], n)
    
    d <- cnt / read_lengths * 100
    d[!is.finite(d)] <- NA_real_
    out$density[[k]] <- d
    
    has <- cnt > 0
    out$means[[k]] <- sm[has] / (255 * cnt[has])        # per-read mean(values / 255)
    out$sum[k] <- sum(sm)
    out$n[k]   <- sum(cnt)
    if (KEEP_ML_VALUES && length(i) > 0L) {
      out$values[[k]] <- ml_all[sequence(ncalls[i], from = gstart[i])]
    }
  }
  out
}

# ---- main function ----------------------------------------------------------

process_file_report <- function(filename, qc_obj) {
  
  sourceCpp("./C/weinr.cpp")
  
  CHUNK_SIZE <- 10000L
  MAX_WIDTH <- 2000000L
  BIN_WIDTH <- 5
  BINS <- 20
  
  per_read_avg_q_score <- c()
  chunk_qpp_subtotals <- vector(mode = "list")
  qpp_frequency_table <- matrix(0, nrow = BINS, ncol = MAX_WIDTH)
  
  if (grepl("\\.bam$", filename, ignore.case = TRUE)) {
    
    bam_stream <- BamFile(filename, yieldSize = CHUNK_SIZE)
    open(bam_stream)
    
    param <- ScanBamParam(
      what = c("seq", "qual"),
      tag = c("pt", "MM", "ML")
    )
    
    polyA_chunks   <- list()
    density_chunks <- list()
    means_chunks   <- list()
    values_chunks  <- list()
    tot_sum <- setNames(numeric(n_mod), mod_names)
    tot_n   <- setNames(numeric(n_mod), mod_names)
    
  } else {
    
    fastq_stream <- FastqStreamer(filename, n = CHUNK_SIZE)
    
  }
  
  start_job <- Sys.time()
  
  per_read_lengths <- numeric(0)
  per_read_N_content <- numeric(0)
  per_read_GC_content <- numeric(0)
  per_read_avg_q_score <- numeric(0)
  
  qpp_count <- numeric(MAX_WIDTH)
  qpp_sum <- numeric(MAX_WIDTH)
  ci <- 0L
  
  repeat {
    
    if (grepl("\\.bam$", filename, ignore.case = TRUE)) {
      
      bam_chunk <- scanBam(bam_stream, param = param)[[1]]
      
      if (length(bam_chunk$seq) == 0) {
        break
      }
      
      ci <- ci + 1L
      
      chunk_qsds <- QualityScaledDNAStringSet(
        DNAStringSet(bam_chunk$seq),
        PhredQuality(bam_chunk$qual)
      )
      
      # BAM-SPECIFIC PROCESSING
      
      read_lengths <- Biostrings::width(chunk_qsds)
      
      pt <- as.numeric(tag_first(bam_chunk$tag$pt, NA_real_))
      polyA_chunks[[ci]] <- pt[is.finite(pt) & pt >= 0]
      
      res <- parse_chunk(
        bam_chunk$tag$MM,
        bam_chunk$tag$ML,
        read_lengths
      )
      
      density_chunks[[ci]] <- res$density
      means_chunks[[ci]]   <- res$means
      values_chunks[[ci]]  <- res$values
      
      tot_sum <- tot_sum + res$sum
      tot_n   <- tot_n + res$n
      
    } else {
      
      fastq_chunk <- yield(fastq_stream)
      
      if (length(fastq_chunk) == 0) {
        break
      }
      
      chunk_qsds <- as(fastq_chunk, "QualityScaledDNAStringSet")
      
    }
    
    # AGGREGATE DNA METRICS
    
    chunk_dna_ss <- Biostrings::DNAStringSet(chunk_qsds)
    
    chunk_per_read_lengths <- Biostrings::width(chunk_dna_ss)
    per_read_lengths <- c(per_read_lengths, chunk_per_read_lengths)
    
    chunk_per_read_base_frequency <- Biostrings::alphabetFrequency(chunk_dna_ss)
    
    chunk_per_read_N_content <- chunk_per_read_base_frequency[, 'N']
    per_read_N_content <- c(per_read_N_content, chunk_per_read_N_content)
    
    chunk_per_read_G_content <- chunk_per_read_base_frequency[, 'G']
    chunk_per_read_C_content <- chunk_per_read_base_frequency[, 'C']
    chunk_per_read_GC_content <- (chunk_per_read_G_content + chunk_per_read_C_content) / chunk_per_read_lengths
    per_read_GC_content <- c(per_read_GC_content, chunk_per_read_GC_content)
    
    # AGGREGATE QUALITY METRICS
    
    chunk_quality_ss <- as.character(Biostrings::quality(chunk_qsds))
    chunk_error_p <- sanger_to_probability(chunk_quality_ss)
    
    chunk_per_read_avg_q_score <- calculate_per_read_avg_q_score(chunk_error_p)
    per_read_avg_q_score <- c(per_read_avg_q_score, chunk_per_read_avg_q_score)
    
    chunk_qpp_subtotals <- tabulate_qpp_mean(chunk_error_p, MAX_WIDTH)
    qpp_count <- qpp_count + chunk_qpp_subtotals[[1]]
    qpp_sum <- qpp_sum + chunk_qpp_subtotals[[2]]
    
    qpp_frequency_table <- tabulate_qpp_median(qpp_frequency_table, chunk_quality_ss, BIN_WIDTH)
    
  }
  
  # CLOSE STREAM
  
  if (grepl("\\.bam$", filename, ignore.case = TRUE)) {
    
    close(bam_stream)
    
    # COMBINE BAM-SPECIFIC METRICS
    
    polyA_lengths <- unlist(polyA_chunks, use.names = FALSE)
    
    collect <- function(chunks) {
      setNames(lapply(mod_names, function(m)
        unlist(lapply(chunks, `[[`, m), use.names = FALSE)), mod_names)
    }
    
    modification_density <- collect(density_chunks)
    read_means <- collect(means_chunks)
    type_values <- if (KEEP_ML_VALUES) collect(values_chunks) else NULL
    
    average_ML_file <- as.list(ifelse(tot_n > 0, tot_sum / (255 * tot_n), NA_real_))
    
  } else {
    
    close(fastq_stream)
    
  }
  
  # CALCULATE FILE METRICS
  
  yield <- sum(per_read_lengths)
  read_count <- length(per_read_lengths)
  N_count <- sum(per_read_N_content)
  
  if (read_count == 0) {
    stop("No reads found in file: ", filename)
  }
  
  decreasinglengths <- sort(per_read_lengths, decreasing = TRUE)
  tmp <- cumsum(as.numeric(decreasinglengths))
  N50 <- decreasinglengths[which(tmp >= yield/2)[1]]
  N90 <- decreasinglengths[which(tmp >= yield*0.9)[1]]
  
  avg_Q_score <- -10 * log10(sum(qpp_sum) / sum(qpp_count))
  qpp_mean <- calculate_qpp_mean(qpp_count, qpp_sum, MAX_WIDTH)
  qpp_median <- calculate_qpp_median(qpp_frequency_table, MAX_WIDTH, BINS, BIN_WIDTH)
  
  end_job <- Sys.time()
  print(end_job - start_job)
  
  # STORE METRICS IN qc_obj
  
  MAX_LENGTH <- max(per_read_lengths)
  mean <- qpp_mean[1:MAX_LENGTH]
  median <- qpp_median[1:MAX_LENGTH]
  
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
  
  if (grepl("\\.bam$", filename, ignore.case = TRUE)) {
    
    qc_obj@metrics[[filename]]$polyALengths <- polyA_lengths
    qc_obj@metrics[[filename]]$modDensity   <- modification_density
    qc_obj@metrics[[filename]]$modReadMeans <- read_means
    qc_obj@metrics[[filename]]$typeValues   <- type_values
    qc_obj@metrics[[filename]]$avgML        <- average_ML_file
    
  }
  
  summary_df <- data.frame(
    File = filename,
    `Read Count` = read_count,
    `Yield (bp)` = yield,
    N50 = N50,
    N90 = N90,
    `Average Per Read Q Score` = avg_Q_score,
    `N Count` = N_count,
    stringsAsFactors = FALSE
  )
  
  qc_obj@summary_metrics <- rbind(qc_obj@summary_metrics, summary_df)
  
  qc_obj
}