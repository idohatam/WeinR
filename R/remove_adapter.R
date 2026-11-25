

remove_adapter <- function(chunk_qsds, 
                           single_qc,
                           adapter_summary,
                           adapter_sequence,
                           max_mismatched_end = 3,
                           min_overlapped_end = 20,
                           min_internal_distance = 100,
                           min_fragment_length = 100) {
  
  if (!is.character(adapter_sequence) || length(adapter_sequence) != 1L || is.na(adapter_sequence) || !nzchar(adapter_sequence))
    stop("RemoveAdapter(): adapterSeq must be non-empty character(1).", call. = FALSE)

  to_int1 <- function(x, nm, min = 0L) {
    if (is.null(x) || !is.numeric(x) || length(x) != 1L || is.na(x) || x < min || x %% 1 != 0)
      stop("RemoveAdapter(): ", nm, " must be a single integer >= ", min, ".", call. = FALSE)
    as.integer(x)
  }
  max_mismatched_end <- to_int1(max_mismatched_end, "MaxMismatchEnd", 0L)
  min_overlapped_end  <- to_int1(min_overlapped_end, "MinOverlapEnd",  0L)
  min_internal_distance <- to_int1(min_internal_distance, "MinInternalDistance", 0L)
  min_fragment_length   <- to_int1(min_fragment_length, "MinFragmentLength",   1L)

  n_reads <- length(chunk_qsds)
  widths <- Biostrings::width(chunk_qsds)
  key <- basename(qc_obj@files[[1]])

  # vectorized matching (both orientations)
  adapter <- Biostrings::DNAString(adapter_sequence)
  adapter_rc <- Biostrings::reverseComplement(adapter)

  m_fwd <- Biostrings::vmatchPattern(adapter, chunk_qsds, max.mismatch = max_mismatched_end)
  m_rev <- Biostrings::vmatchPattern(adapter_rc, chunk_qsds, max.mismatch = max_mismatched_end)

  # get first start and last end quickly per read (NA if none)
  first_start <- integer(n_reads); first_start[] <- NA_integer_
  first_end   <- integer(n_reads); first_end[]   <- NA_integer_
  last_start  <- integer(n_reads); last_start[]  <- NA_integer_
  last_end    <- integer(n_reads); last_end[]    <- NA_integer_

  has_any <- logical(n_reads)
  
  for (i in seq_len(n_reads)) {
    h <- c(m_fwd[[i]], m_rev[[i]])
    if (length(h) == 0L) next

    # stable ordering without sort() dispatch
    h <- h[order(IRanges::start(h), IRanges::end(h))]

    has_any[i] <- TRUE
    first_start[i] <- IRanges::start(h)[1]
    first_end[i]   <- IRanges::end(h)[1]
    last_start[i]  <- IRanges::start(h)[length(h)]
    last_end[i]    <- IRanges::end(h)[length(h)]
  }
  
  # vectorized END trimming
  trim_start <- rep.int(1L, n_reads)
  trim_end   <- widths

  has_5prime <- !is.na(first_start) & first_start <= min_overlapped_end
  trim_start[has_5prime] <- pmin(first_end[has_5prime] + 1L, widths[has_5prime] + 1L)

  has_3prime <- !is.na(last_end) & (widths - last_end) <= min_overlapped_end
  trim_end[has_3prime] <- pmax(last_start[has_3prime] - 1L, 1L)

  keep_endtrim <- trim_end >= trim_start
  if (!any(keep_endtrim)) {
    warning("RemoveAdapter(): all reads became invalid after end trimming for: ", key)
    qc_obj@metadata[[key]]$adapter_summary <- list(
      reads_before = n_reads, reads_after = 0L, output_paths = character(0)
    )
    return(qc_obj)
  }
  
  # Use subseq() for stability; it works on QualityScaledDNAStringSet too
  endtrimmed <- Biostrings::subseq(chunk_qsds[keep_endtrim],
                                   start = trim_start[keep_endtrim],
                                   end   = trim_end[keep_endtrim])

  # Enforce min length after end trim
  endtrimmed <- endtrimmed[Biostrings::width(endtrimmed) >= min_fragment_length]
  
  # internal splitting: loop only over reads with >=2 hits
  # We do strict internal search (exact) like your old code, but only for candidates.
  chimeric_frags <- chunk_qsds[0]   # safe empty
  
  if (min_internal_distance > 0L) {

    # candidate reads: those with >=2 end-match hits in either direction (cheap heuristic)
    cand <- which(vapply(seq_len(n_reads), function(i) {
      length(c(m_fwd[[i]], m_rev[[i]])) >= 2L
    }, logical(1)))

    if (length(cand) > 0L) {

      qual_set <- Biostrings::quality(chunk_qsds)

      for (i in cand) {
        seq_i  <- chunk_qsds[[i]]
        qual_i <- qual_set[[i]]
        L <- widths[i]

        # strict internal exact hits (both orientations)
        h_mid <- c(
          Biostrings::matchPattern(adapter,    seq_i, max.mismatch = 0),
          Biostrings::matchPattern(adapter_rc, seq_i, max.mismatch = 0)
        )
        if (length(h_mid) < 1L) next
        h_mid <- h_mid[order(IRanges::start(h_mid), IRanges::end(h_mid))]

        valid_internal <- which(
          IRanges::start(h_mid) >= min_internal_distance &
            IRanges::end(h_mid) <= (L - min_internal_distance)
        )
        if (length(valid_internal) == 0L) next

        read_name <- names(chunk_qsds)[i]
        if (is.null(read_name) || is.na(read_name) || !nzchar(read_name))
          read_name <- paste0("read", i)

        if (TRUE) message("Splitting chimeric read ", read_name)

        split_starts <- c(1L, IRanges::end(h_mid)[valid_internal] + 1L)
        split_ends   <- c(IRanges::start(h_mid)[valid_internal] - 1L, L)

        for (frag_idx in seq_along(split_starts)) {
          s <- split_starts[frag_idx]
          e <- split_ends[frag_idx]
          if (e < s) next
          if ((e - s + 1L) < min_fragment_length) next

          frag_seq_raw  <- Biostrings::subseq(seq_i,  start = s, end = e)
          frag_qual_raw <- Biostrings::subseq(qual_i, start = s, end = e)

          frag_seq <- Biostrings::DNAStringSet(frag_seq_raw)
          frag_qual <- Biostrings::PhredQuality(Biostrings::BStringSet(as.character(frag_qual_raw)))

          nm <- paste0(read_name, "_part", frag_idx)
          names(frag_seq)  <- nm
          names(frag_qual) <- nm

          chimeric_frags <- c(chimeric_frags, Biostrings::QualityScaledDNAStringSet(frag_seq, frag_qual))
        }
      }
    }
  }

  # Remove originals that were split, then add fragments
  if (length(chimeric_frags) > 0L) {
    originals <- unique(sub("_part.*", "", names(chimeric_frags)))
    drop_idx <- which(sub("_part.*", "", names(endtrimmed)) %in% originals)
    if (length(drop_idx) > 0L) endtrimmed <- endtrimmed[-drop_idx]
  }

  final_reads <- c(endtrimmed, chimeric_frags)
  # final_reads <- final_reads[Biostrings::width(final_reads) >= min_fragment_length]

  n_after <- length(final_reads)
  adapter_summary$reads_before <- adapter_summary$reads_before + n_reads
  adapter_summary$reads_after <- adapter_summary$reads_after + n_after
    
  final_reads
}