

validate_input <- function(files, minimum_reads) {
  
  if (!is.numeric(minimum_reads) || length(minimum_reads) != 1L || is.na(minimum_reads) ||
      minimum_reads < 1 || minimum_reads %% 1 != 0) {
    stop("minimum_reads must be a single positive integer (>= 1).",
         call. = FALSE
         )
  }
  minimum_reads <- as.integer(minimum_reads)
  
  if (length(files) == 1L && dir.exists(files)) {
    input_dir <- normalizePath(files, mustWork = TRUE)
    files <- list.files(
      path = input_dir,
      pattern = "\\.(fastq|fastq\\.gz|bam)$",
      full.names = TRUE
    )
    if (length(files) == 0L) {
      stop(
        "No .fastq, .fastq.gz, or .bam files found in: ",
        input_dir,
        call. = FALSE
      )
    }
  }
    
  for(file in files) {
    file_type <- CheckFile(file)
    print(file_type)
    if (file_type %in% c("fastq", "fastq.gz")) {
      tryCatch(FastqCheck(file, n_records = 100L),
               error = function(e) stop("FASTQ check failed (corrupted/malformed): ", conditionMessage(e), call. = FALSE))
    } else if (file_type == "bam") {
      bf <- Rsamtools::BamFile(file, yieldSize = 1000)
      tryCatch(
        bam <- Rsamtools::scanBam(bf)[[1]],
        error = function(e) {
          stop(conditionMessage(e), call. = FALSE)
        }
      )
      if (is.null(bam$seq) || is.null(bam$qual)) {
        stop("BAM missing required fields 'seq' and/or 'qual'.", call. = FALSE)
      }
      if (length(bam$seq) != length(bam$qual)) {
        stop("BAM 'seq' and 'qual' lengths do not match.", call. = FALSE)
      }
    }
  }
  qc_obj <- .init_qc_object(files)
  
}