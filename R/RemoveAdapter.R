#' Trim adapters and split chimeric reads (internal)
#'
#' Internal helper that searches for an adapter sequence within reads from a
#' single-file QC object, trims adapter occurrences near read ends, and optionally
#' splits reads with internal adapter hits into fragments (chimeric read handling).
#' Trimmed reads are written to disk and a summary is stored in
#' \code{qc_obj@metadata[[basename(original_file)]]}.
#'
#' @param qc_obj A \code{LongReadQC} (or compatible) object containing exactly one
#'   file in \code{@files}. This function operates on a single file at a time.
#' @param adapterSeq Character(1). Adapter sequence (DNA alphabet) to search for.
#' @param MaxMismatchEnd Integer(1). Maximum mismatches allowed when matching the
#'   adapter sequence.
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
#' @param OutFile Character(1) or \code{NULL}. Optional output filename stem used
#'   for naming outputs. If \code{NULL}, naming is based on the original input file.
#' @param verbose Logical(1). If \code{TRUE}, prints progress messages.
#' @param WriteIntermediate Logical(1). If \code{TRUE}, writes an intermediate FASTQ
#'   file and stores the path in \code{attr(qc_obj, "._tmp_fastq")} for chaining
#'   steps in downstream workflows.
#'
#' @return The updated \code{qc_obj}. An adapter trimming summary is stored in
#'   \code{qc_obj@metadata[[basename(original_file)]]$adapter_summary}.
#'
#' @details
#' Adapter detection uses \code{Biostrings::vmatchPattern()} with
#' \code{max.mismatch = MaxMismatchEnd}. Adapter hits within \code{MinOverlapEnd}
#' of the 5' or 3' ends are treated as end adapters and trimmed.
#'
#' Reads with two or more adapter hits may contain internal adapters. Internal hits
#' that are at least \code{MinInternalDistance} from both ends are used to split a
#' read into fragments (chimera handling). Fragments shorter than
#' \code{MinFragmentLength} are discarded. If chimeric fragments are produced, the
#' corresponding original reads are removed from the end-trimmed set to avoid
#' duplicate representation.
#'
#' The summary is recorded under the original input file (the single file stored in
#' \code{qc_obj@files[[1]]}), keyed by its basename.
#'
#' @examples
#' \dontrun{
#' qc_obj2 <- RemoveAdapter(
#'   qc_obj,
#'   adapterSeq = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
#'   MaxMismatchEnd = 3,
#'   MinOverlapEnd = 20,
#'   MinInternalDistance = 100,
#'   MinFragmentLength = 200,
#'   OutDir = tempdir(),
#'   OutFileType = "fastq",
#'   verbose = TRUE,
#'   WriteIntermediate = TRUE
#' )
#' }
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
                          OutFile = NULL,
                          verbose = TRUE,
                          WriteIntermediate = FALSE) {
  
  ## Input checks
  if (length(qc_obj@files) != 1)
    stop("RemoveAdapter() expects a LongReadQC object with exactly one file.")
  
  attr(qc_obj, "._tmp_fastq") <- NULL
  
  # Choose which file to trim
  fpath <- if (!is.null(FilePath)) FilePath else qc_obj@files[[1]]
  OriginalFPath <- qc_obj@files[[1]]
  fname <- basename(fpath)
  
  if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  if (verbose) message("Adapter trimming: ", fname)
  
  ## Load reads
  reads <- ImportFile(fpath)
  
  if (!inherits(reads, "QualityScaledDNAStringSet"))
    stop("ImportFile() must return a QualityScaledDNAStringSet object.")
  
  adapter <- Biostrings::DNAString(adapterSeq)
  read_widths <- Biostrings::width(reads)
  n_reads <- length(reads)
  
  ## Vectorized adapter hit detection
  matches <- Biostrings::vmatchPattern(
    adapter, reads, max.mismatch = MaxMismatchEnd
  )
  
  ## extract first hit start / last hit end per read
  starts <- vapply(matches, function(h)
    if (length(h)) Biostrings::start(h)[1] else NA_integer_,
    integer(1)
  )
  
  ends <- vapply(matches, function(h)
    if (length(h)) Biostrings::end(h)[length(h)] else NA_integer_,
    integer(1)
  )
  
  
  ## END TRIMMING
  trim_start <- rep(1L, n_reads)
  trim_end   <- read_widths
  
  ## 5' trimming
  has_5prime <- !is.na(starts) & starts <= MinOverlapEnd
  trim_start[has_5prime] <- ends[has_5prime] + 1L
  
  ## 3' trimming
  has_3prime <- !is.na(ends) & (read_widths - ends) <= MinOverlapEnd
  trim_end[has_3prime] <- pmax(starts[has_3prime] - 1L, 1L)
  
  ## Apply end trimming
  trimmed_reads <- IRanges::narrow(reads, start = trim_start, end = trim_end)
  
  
  ## internal adapter trimming
  chimeric_list <- list()
  
  for (i in seq_len(n_reads)) {
    
    hit_coords <- as.data.frame(matches[[i]])
    if (nrow(hit_coords) < 2) next  
    
    seq_i  <- reads[[i]]
    qual_i <- Biostrings::quality(reads)[[i]]
    seq_len <- Biostrings::nchar(seq_i)
    read_name <- names(reads)[i]
    
    ## filter hits far from ends
    valid_hits <- subset(
      hit_coords,
      start >= MinInternalDistance &
        end <= (seq_len - MinInternalDistance)
    )
    if (nrow(valid_hits) == 0) next
    
    ## define split intervals
    split_starts <- c(1, valid_hits$end + 1)
    split_ends   <- c(valid_hits$start - 1, seq_len)
    
    frag_lengths <- split_ends - split_starts + 1
    keep_idx <- which(frag_lengths >= MinFragmentLength)
    if (length(keep_idx) == 0) next
    
    ## extract fragment sequences + qualities
    frags <- mapply(function(s, e, idx) {
      frag_seq_raw  <- Biostrings::subseq(seq_i, s, e)
      frag_qual_raw <- Biostrings::subseq(qual_i, s, e)
      
      frag_seq <- Biostrings::DNAStringSet(frag_seq_raw)
      frag_qual <- Biostrings::PhredQuality(
        Biostrings::BStringSet(as.character(frag_qual_raw))
      )
      
      names(frag_seq)  <- paste0(read_name, "_part", idx)
      names(frag_qual) <- paste0(read_name, "_part", idx)
      
      Biostrings::QualityScaledDNAStringSet(frag_seq, frag_qual)
    },
    split_starts[keep_idx],
    split_ends[keep_idx],
    seq_along(keep_idx),
    SIMPLIFY = FALSE)
    
    chimeric_list <- c(chimeric_list, frags)
  }
  
  
  if (length(chimeric_list) > 0) {
    chimeric_frags <- do.call(c, chimeric_list)
  } else {
    chimeric_frags <- Biostrings::QualityScaledDNAStringSet()
  }
  
  
  if (length(chimeric_frags) > 0) {
    
    ## extract original names
    original_names <- unique(sub("_part.*", "", names(chimeric_frags)))
    
    ## find them in trimmed_reads
    to_drop <- which(sub("_part.*", "", names(trimmed_reads)) %in% original_names)
    
    if (length(to_drop) > 0) {
      trimmed_reads <- trimmed_reads[-to_drop]
    }
  }
  
  
  ## Combine final outputs
  all_trimmed <- c(trimmed_reads, chimeric_frags)
  
  ## keep only fragments >= MinFragmentLength
  keep <- Biostrings::width(all_trimmed) >= MinFragmentLength
  all_trimmed <- all_trimmed[keep]
  
  dna <- Biostrings::DNAStringSet(all_trimmed)
  qual <- Biostrings::quality(all_trimmed)
  names(qual) <- names(dna)
  
  all_trimmed <- Biostrings::QualityScaledDNAStringSet(dna, qual)
  
  
  ## Write output
  base_name <- if (!is.null(OutFile))
    RemoveExt(OutFile) else RemoveExt(OriginalFPath)
  
  out_paths <- WriteReadOutputs(
    reads = all_trimmed,
    base_name = paste0(base_name, "_adaptertrimmed"),
    OutDir = OutDir,
    OutFileType = OutFileType
  )
  
  
  ## Update metadata
  qc_obj@metadata[[basename(OriginalFPath)]]$adapter_summary <- list(
    reads_before = n_reads,
    reads_after  = length(all_trimmed),
    output_paths = out_paths
  )
  
  if (verbose) {
    message("Adapter trimming done for ", fname,
            ": kept ", length(all_trimmed), " reads (",
            length(chimeric_frags), " chimeric fragments).")
  }
  
  tmp_fastq <- NULL
  
  if (isTRUE(WriteIntermediate) && length(all_trimmed) > 0L) {
    
    stem <- RemoveExt(basename(OriginalFPath))
    tmp_base <- paste0(stem, "_adapter_tmp_", Sys.getpid(), "_", sample.int(1e9, 1))
    
    tmp_dir <- file.path(tempdir(), "longR_intermediate")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    
    tmp_paths <- WriteReadOutputs(
      reads = all_trimmed,
      base_name = tmp_base,
      OutDir = tmp_dir,
      OutFileType = "fastq",
      Verbose = FALSE
    )
    
    tmp_fastq <- tmp_paths[[1]]
  }
  
  
  attr(qc_obj, "._tmp_fastq") <- tmp_fastq
  
  base::return(qc_obj)
}

