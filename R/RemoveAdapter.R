#' Trim adapters and split chimeric reads (internal)
#'
#' Internal helper that searches for an adapter sequence within reads from a
#' single-file QC object, trims adapter occurrences near read ends, and optionally
#' splits reads with internal adapter hits into fragments (chimeric read handling).
#' A summary is stored in \code{qc_obj@metadata[[basename(original_file)]]$adapter_summary}.
#' Optionally, an intermediate FASTQ can be written for chaining and stored in
#' \code{attr(qc_obj, "._tmp_fastq")}.
#'
#' @param qc_obj A \code{LongReadQC} (or compatible) object containing exactly one
#'   file in \code{@files}. This function operates on a single file at a time.
#' @param adapterSeq Character(1). Adapter sequence (DNA alphabet) to search for.
#' @param MaxMismatchEnd Integer(1). Maximum mismatches allowed when matching the
#'   adapter sequence during end-adapter detection.
#' @param MinOverlapEnd Integer(1). Maximum distance from either end (5' or 3')
#'   within which an adapter hit is treated as an end adapter and trimmed.
#' @param MinInternalDistance Integer(1). Minimum distance from either end required
#'   for an adapter hit to be considered internal (used for chimera splitting).
#' @param MinFragmentLength Integer(1). Minimum fragment length to keep after
#'   trimming/splitting. Fragments shorter than this are discarded.
#' @param FilePath Character(1) or \code{NULL}. Optional file path to trim. If
#'   \code{NULL}, uses the file in \code{qc_obj@files[[1]]}. This supports trimming
#'   of intermediate files while still recording results under the original input file.
#' @param OutDir Character(1). Output directory. Created if it does not exist.
#' @param OutFileType Character vector specifying output format(s). Supported values
#'   are \code{"fastq"}, \code{"fasta"}, and \code{"bam"} (as implemented by
#'   \code{WriteReadOutputs()}). Multiple formats may be requested.
#' @param verbose Logical(1). If \code{TRUE}, prints progress messages.
#' @param WriteIntermediate Logical(1). If \code{TRUE}, writes an intermediate FASTQ
#'   file and stores the path in \code{attr(qc_obj, "._tmp_fastq")} for chaining.
#' @param KeepIntermediates Logical(1). If \code{TRUE}, writes the processed reads to
#'   \code{OutDir} (in the formats requested by \code{OutFileType}) and records the
#'   written paths in the summary. If \code{FALSE}, no final outputs are written.
#' @param OutSuffix Character(1). Suffix appended to the output basename when writing
#'   processed reads (e.g., \code{"adaptertrimmed"}).
#'
#' @return The updated \code{qc_obj}. An adapter trimming summary is stored in
#'   \code{qc_obj@metadata[[basename(original_file)]]$adapter_summary}. If
#'   \code{WriteIntermediate = TRUE}, \code{attr(qc_obj, "._tmp_fastq")} contains the
#'   intermediate FASTQ path (or \code{NULL}).
#'
#' @details
#' Adapter detection uses \code{Biostrings::vmatchPattern()} (both orientations:
#' adapter and reverse complement) with \code{max.mismatch = MaxMismatchEnd}.
#' Adapter hits within \code{MinOverlapEnd} of the 5' or 3' ends are treated as end
#' adapters and trimmed.
#'
#' Reads with internal adapter hits at least \code{MinInternalDistance} from both ends
#' can be split into fragments. Fragments shorter than \code{MinFragmentLength} are
#' discarded. If fragments are produced, corresponding original reads are removed from
#' the end-trimmed set to avoid duplicates.
#'
#' @keywords internal
RemoveAdapter <- function(qc_obj,
                          adapterSeq,
                          MaxMismatchEnd = 3,
                          MinOverlapEnd = 20,
                          MinInternalDistance = 100,
                          MinFragmentLength = 200,
                          FilePath = NULL,
                          OutDir = ".",
                          OutFileType = c("fastq"),
                          verbose = TRUE,
                          WriteIntermediate = FALSE,
                          KeepIntermediates = FALSE,
                          OutSuffix = "adaptertrimmed") {
  
  # ---- checks ----
  if (length(qc_obj@files) != 1L)
    stop("RemoveAdapter() expects a QC object with exactly one file.")
  if (!is.character(adapterSeq) || length(adapterSeq) != 1L || is.na(adapterSeq) || !nzchar(adapterSeq))
    stop("RemoveAdapter(): adapterSeq must be non-empty character(1).")
  
  to_int1 <- function(x, nm, min = 0L) {
    if (is.null(x) || !is.numeric(x) || length(x) != 1L || is.na(x) || x < min || x %% 1 != 0)
      stop("RemoveAdapter(): ", nm, " must be a single integer >= ", min, ".")
    as.integer(x)
  }
  MaxMismatchEnd <- to_int1(MaxMismatchEnd, "MaxMismatchEnd", 0L)
  MinOverlapEnd  <- to_int1(MinOverlapEnd,  "MinOverlapEnd",  0L)
  MinInternalDistance <- to_int1(MinInternalDistance, "MinInternalDistance", 0L)
  MinFragmentLength   <- to_int1(MinFragmentLength,   "MinFragmentLength",   1L)
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)
  attr(qc_obj, "._tmp_fastq") <- NULL
  
  # ---- choose file (supports chaining) ----
  fpath <- if (!is.null(FilePath)) FilePath else qc_obj@files[[1]]
  OriginalFPath <- qc_obj@files[[1]]
  key <- basename(OriginalFPath)
  fname <- basename(fpath)
  if (isTRUE(verbose)) message("Adapter trimming: ", key)
  
  reads <- ImportFile(fpath)
  if (!inherits(reads, "QualityScaledDNAStringSet"))
    stop("RemoveAdapter(): ImportFile() must return a QualityScaledDNAStringSet.")
  
  n_reads <- length(reads)
  if (n_reads == 0L) {
    warning("RemoveAdapter(): input contains 0 reads for: ", key)
    qc_obj@metadata[[key]]$adapter_summary <- list(
      reads_before = 0L, reads_after = 0L, output_paths = character(0)
    )
    return(qc_obj)
  }
  
  widths <- Biostrings::width(reads)
  
  # ---- vectorized matching (both orientations) ----
  adapter <- Biostrings::DNAString(adapterSeq)
  adapter_rc <- Biostrings::reverseComplement(adapter)
  
  m_fwd <- Biostrings::vmatchPattern(adapter,    reads, max.mismatch = MaxMismatchEnd)
  m_rev <- Biostrings::vmatchPattern(adapter_rc, reads, max.mismatch = MaxMismatchEnd)
  
  # helper: get first start and last end quickly per read (NA if none)
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
  
  # ---- vectorized END trimming ----
  trim_start <- rep.int(1L, n_reads)
  trim_end   <- widths
  
  has_5prime <- !is.na(first_start) & first_start <= MinOverlapEnd
  trim_start[has_5prime] <- pmin(first_end[has_5prime] + 1L, widths[has_5prime] + 1L)
  
  has_3prime <- !is.na(last_end) & (widths - last_end) <= MinOverlapEnd
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
  endtrimmed <- Biostrings::subseq(reads[keep_endtrim],
                                   start = trim_start[keep_endtrim],
                                   end   = trim_end[keep_endtrim])
  
  # Enforce min length after end trim
  endtrimmed <- endtrimmed[Biostrings::width(endtrimmed) >= MinFragmentLength]
  
  # ---- internal splitting: loop only over reads with >=2 hits ----
  # We do strict internal search (exact) like your old code, but only for candidates.
  chimeric_frags <- reads[0]   # safe empty
  
  if (MinInternalDistance > 0L) {
    
    # candidate reads: those with >=2 end-match hits in either direction (cheap heuristic)
    cand <- which(vapply(seq_len(n_reads), function(i) {
      length(c(m_fwd[[i]], m_rev[[i]])) >= 2L
    }, logical(1)))
    
    if (length(cand) > 0L) {
      
      qual_set <- Biostrings::quality(reads)
      
      for (i in cand) {
        seq_i  <- reads[[i]]
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
          IRanges::start(h_mid) >= MinInternalDistance &
            IRanges::end(h_mid) <= (L - MinInternalDistance)
        )
        if (length(valid_internal) == 0L) next
        
        read_name <- names(reads)[i]
        if (is.null(read_name) || is.na(read_name) || !nzchar(read_name))
          read_name <- paste0("read", i)
        
        if (isTRUE(verbose)) message("Splitting chimeric read ", read_name)
        
        split_starts <- c(1L, IRanges::end(h_mid)[valid_internal] + 1L)
        split_ends   <- c(IRanges::start(h_mid)[valid_internal] - 1L, L)
        
        for (frag_idx in seq_along(split_starts)) {
          s <- split_starts[frag_idx]
          e <- split_ends[frag_idx]
          if (e < s) next
          if ((e - s + 1L) < MinFragmentLength) next
          
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
  final_reads <- final_reads[Biostrings::width(final_reads) >= MinFragmentLength]
  
  n_after <- length(final_reads)
  drop_frac <- (n_reads - n_after) / n_reads
  
  if (n_after == 0L) {
    warning("RemoveAdapter(): no reads remained after adapter trimming for: ", key)
    qc_obj@metadata[[key]]$adapter_summary <- list(
      reads_before = n_reads, reads_after = 0L, output_paths = character(0)
    )
    return(qc_obj)
  }
  
  if (drop_frac >= 0.50) {
    warning("RemoveAdapter(): trimming removed ",
            round(drop_frac * 100, 1), "% of reads for: ", key,
            " (", n_after, "/", n_reads, " kept).")
  }
 
  
  out_paths <- character(0)
  if (isTRUE(KeepIntermediates)) {
    base_name <- paste0(RemoveExt(OriginalFPath), "_", OutSuffix)
    out_paths <- WriteReadOutputs(
      reads = final_reads,
      BaseName = base_name,
      OutDir = OutDir,
      OutFileType = OutFileType
    )
  }
  qc_obj@metadata[[key]]$adapter_summary <- list(
    reads_before = n_reads,
    reads_after  = n_after,
    output_paths = out_paths
  )
  
  if (isTRUE(verbose)) message("Adapter trimming done for ", key, ": kept ", n_after, " reads.")
  
  tmp_fastq <- NULL
  if (isTRUE(WriteIntermediate)) {
    stem <- RemoveExt(basename(OriginalFPath))
    tmp_base <- paste0(stem, "_adapter_tmp_", Sys.getpid(), "_", sample.int(1e9, 1))
    tmp_dir <- file.path(tempdir(), "longR_intermediate")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    tmp_paths <- WriteReadOutputs(final_reads, tmp_base, tmp_dir, "fastq", Verbose = FALSE)
    tmp_fastq <- tmp_paths[[1]]
  }
  
  attr(qc_obj, "._tmp_fastq") <- tmp_fastq
  qc_obj
}
