make_qc_onefile <- function(fpath, readLengths, meanprQscore, Ncount) {
  key <- basename(fpath)
  
  qc <- .init_qc_object(files = fpath)
  
  # Add what downstream expects
  qc@metrics <- list()
  qc@metrics[[key]] <- list(
    readLengths  = readLengths,
    meanprQscore = meanprQscore,
    Ncount       = Ncount
  )
  
  # Make sure per-file metadata key exists
  qc@metadata[[key]] <- list()
  
  # Make summary_metrics$file match (ProcessReads subsets on it)
  qc@summary_metrics <- data.frame(
    file = key,
    yield = 0, N50 = 0, N90 = 0, avgQscore = 0, `N count` = 0,
    stringsAsFactors = FALSE
  )
  
  qc
}

fake_reads <- function(n = 4, width = 20) {
  testthat::skip_if_not_installed("Biostrings")
  seqs  <- rep(paste(rep("A", width), collapse = ""), n)
  quals <- rep(paste(rep("F", width), collapse = ""), n)
  
  Biostrings::QualityScaledDNAStringSet(
    Biostrings::DNAStringSet(seqs),
    Biostrings::PhredQuality(Biostrings::BStringSet(quals))
  )
}
