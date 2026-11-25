

trim_long <- function(chunk_qsds, single_qc, trim_summary, start_index, end_index, start = NULL, end = NULL) {
  
  filename <- single_qc@files[[1]]
  key <- basename(filename)
  
  if (length(qc_obj@files) != 1L) {
    stop("TrimLong() expects a QC object with exactly one file in qc_obj@files.", call. = FALSE)
  }
  validate_endpoint <- function(position, name) {
    if (is.null(position)) return(invisible(TRUE))
    if (!is.numeric(position) || length(position) != 1L || is.na(position) || position < 0 || position %% 1 != 0) {
      stop("TrimLong(): ", name, " must be NULL or a single non-negative integer.", call. = FALSE)
    }
    invisible(TRUE)
  }
  validate_endpoint(start, "Start")
  validate_endpoint(end, "End")
  
  n_reads <- length(chunk_qsds)
  n <- Biostrings::width(chunk_qsds)
  Start_i <- if (!is.null(start)) as.integer(start) else NULL
  End_i   <- if (!is.null(end))   as.integer(end)   else NULL
  trim_start <- rep(if (!is.null(Start_i)) Start_i + 1L else 1L, n_reads)
  trim_end   <- if (!is.null(End_i)) n - End_i else n
  keep <- trim_end >= trim_start
  n_after <- sum(keep)
  if (n_after == 0L) {
    warning("TrimLong(): no reads remain after trimming for: ", key)
  }
  
  processed_chunk <- IRanges::narrow(
    chunk_qsds[keep],
    start = trim_start[keep],
    end   = trim_end[keep]
  )
  
  trim_summary$reads_before <- trim_summary$reads_before + n_reads
  trim_summary$reads_after <- trim_summary$reads_after + n_after
  
  processed_chunk
}