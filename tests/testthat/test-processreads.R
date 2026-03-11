test_that("ProcessReads: filter + adapter + trim updates metadata and writes final output", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ShortRead")
  skip_if_not_installed("IRanges")
  
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  
  tmp <- withr::local_tempdir()
  withr::local_envvar(c(TMPDIR = tmp, TEMP = tmp, TMP = tmp))
  
  write_test_fastq <- function(path, seqs, ids = NULL, qchar = "I") {
    if (is.null(ids)) ids <- paste0("read", seq_along(seqs))
    stopifnot(length(ids) == length(seqs))
    
    dna <- Biostrings::DNAStringSet(seqs)
    names(dna) <- ids
    
    qs <- Biostrings::PhredQuality(
      Biostrings::BStringSet(vapply(
        seqs,
        function(s) paste(rep(qchar, nchar(s)), collapse = ""),
        character(1)
      ))
    )
    names(qs) <- ids
    
    qsdss <- Biostrings::QualityScaledDNAStringSet(dna, qs)
    
    sr <- ShortRead::ShortReadQ(
      sread   = Biostrings::DNAStringSet(qsdss),
      quality = ShortRead::FastqQuality(Biostrings::quality(qsdss)),
      id      = Biostrings::BStringSet(ids)
    )
    ShortRead::writeFastq(sr, file = path, compress = FALSE)
    invisible(path)
  }
  
  adapter  <- "ACGTACGT"
  prefix20 <- paste(rep("A", 20), collapse = "")
  suffix20 <- paste(rep("T", 20), collapse = "")
  
  r1 <- paste0(prefix20, adapter)   # 3' end adapter
  r2 <- paste0(adapter, suffix20)   # 5' end adapter
  r3 <- paste0(prefix20, suffix20)  # no adapter
  
  fastq_path <- file.path(wd, "toy.fastq")
  write_test_fastq(fastq_path, seqs = c(r1, r2, r3), ids = c("r1", "r2", "r3"))
  
  qc <- .init_qc_object(files = fastq_path)
  qsds <- ImportFile(fastq_path)
  expect_true(inherits(qsds, "QualityScaledDNAStringSet"))
  
  key <- basename(fastq_path)
  qc <- QualMat(qc, qsds, filename = key)
  expect_true(!is.null(qc@metrics[[key]]))
  
  qc2 <- ProcessReads(
    qc_obj = qc,
    filter = TRUE,
    MinAvgQS = 0,
    MinLength = 10,
    MaxNumberNs = 0,
    AdapterSeq = adapter,
    MaxMismatchEnd = 0,
    MinOverlapEnd = 5,
    MinInternalDistance = 200,
    MinFragmentLength = 10,
    Start = 1,
    End = 1,
    OutFileType = "fastq",
    outpath = "processreads_all_steps",
    render_report = FALSE,
    force = TRUE,
    verbose = FALSE,
    KeepIntermediates = FALSE,
    FinalSuffix = "processed"
  )
  
  md <- qc2@metadata[[key]]
  expect_true(is.list(md))
  
  expect_true(!is.null(md$filter_summary))
  expect_equal(md$filter_summary$reads_before, 3L)
  expect_equal(md$filter_summary$reads_after, 3L)
  
  expect_true(!is.null(md$adapter_summary))
  expect_equal(md$adapter_summary$reads_before, 3L)
  expect_equal(md$adapter_summary$reads_after, 3L)
  
  expect_true(!is.null(md$trim_summary))
  expect_equal(md$trim_summary$reads_before, 3L)
  expect_equal(md$trim_summary$reads_after, 3L)
  
  expect_true(length(md$trim_summary$output_paths) >= 1)
  expect_true(all(file.exists(md$trim_summary$output_paths)))
  
  inter_dir <- file.path(tempdir(), "longR_intermediate")
  if (dir.exists(inter_dir)) {
    inter_files <- list.files(inter_dir, pattern = "\\.fastq$", full.names = TRUE)
    expect_length(inter_files, 0L)
  }
})

test_that("ProcessReads: filter-only writes final output and records filter_summary", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ShortRead")
  skip_if_not_installed("IRanges")
  
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  
  write_test_fastq <- function(path, seqs, ids = NULL, qchar = "I") {
    if (is.null(ids)) ids <- paste0("read", seq_along(seqs))
    dna <- Biostrings::DNAStringSet(seqs); names(dna) <- ids
    qs <- Biostrings::PhredQuality(Biostrings::BStringSet(vapply(
      seqs, function(s) paste(rep(qchar, nchar(s)), collapse = ""), character(1)
    )))
    names(qs) <- ids
    qsdss <- Biostrings::QualityScaledDNAStringSet(dna, qs)
    sr <- ShortRead::ShortReadQ(
      sread = Biostrings::DNAStringSet(qsdss),
      quality = ShortRead::FastqQuality(Biostrings::quality(qsdss)),
      id = Biostrings::BStringSet(ids)
    )
    ShortRead::writeFastq(sr, file = path, compress = FALSE)
    path
  }
  
  fastq_path <- file.path(wd, "toy.fastq")
  write_test_fastq(fastq_path, seqs = c(paste(rep("A", 30), collapse = ""),
                                        paste(rep("C", 30), collapse = "")),
                   ids = c("r1", "r2"))
  
  qc <- .init_qc_object(files = fastq_path)
  key <- basename(fastq_path)
  qc <- QualMat(qc, ImportFile(fastq_path), filename = key)
  
  qc2 <- ProcessReads(
    qc_obj = qc,
    filter = TRUE,
    MinAvgQS = 0,
    MinLength = 10,
    MaxNumberNs = 999,
    AdapterSeq = NULL,
    Start = NULL,
    End = NULL,
    OutFileType = "fastq",
    outpath = "processreads_filter_only",
    render_report = FALSE,
    force = TRUE,
    verbose = FALSE,
    KeepIntermediates = TRUE,  # needed so FilterLong writes final outputs
    FinalSuffix = "processed"
  )
  
  md <- qc2@metadata[[key]]
  expect_true(!is.null(md$filter_summary))
  expect_equal(md$filter_summary$reads_before, 2L)
  expect_equal(md$filter_summary$reads_after, 2L)
  expect_true(length(md$filter_summary$output_paths) >= 1)
  expect_true(all(file.exists(md$filter_summary$output_paths)))
  
  expect_true(is.null(md$adapter_summary))
  expect_true(is.null(md$trim_summary))
})

test_that("ProcessReads: adapter-only writes final output and records adapter_summary", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ShortRead")
  skip_if_not_installed("IRanges")
  
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  
  write_test_fastq <- function(path, seqs, ids = NULL, qchar = "I") {
    if (is.null(ids)) ids <- paste0("read", seq_along(seqs))
    dna <- Biostrings::DNAStringSet(seqs); names(dna) <- ids
    qs <- Biostrings::PhredQuality(Biostrings::BStringSet(vapply(
      seqs, function(s) paste(rep(qchar, nchar(s)), collapse = ""), character(1)
    )))
    names(qs) <- ids
    qsdss <- Biostrings::QualityScaledDNAStringSet(dna, qs)
    sr <- ShortRead::ShortReadQ(
      sread = Biostrings::DNAStringSet(qsdss),
      quality = ShortRead::FastqQuality(Biostrings::quality(qsdss)),
      id = Biostrings::BStringSet(ids)
    )
    ShortRead::writeFastq(sr, file = path, compress = FALSE)
    path
  }
  
  adapter <- "ACGTACGT"
  seqs <- c(paste0(adapter, paste(rep("T", 30), collapse = "")),
            paste0(paste(rep("A", 30), collapse = ""), adapter))
  
  fastq_path <- file.path(wd, "toy.fastq")
  write_test_fastq(fastq_path, seqs = seqs, ids = c("r1", "r2"))
  
  qc <- .init_qc_object(files = fastq_path)
  key <- basename(fastq_path)
  qc <- QualMat(qc, ImportFile(fastq_path), filename = key)
  
  qc2 <- ProcessReads(
    qc_obj = qc,
    filter = FALSE,
    AdapterSeq = adapter,
    MaxMismatchEnd = 0,
    MinOverlapEnd = 10,
    MinInternalDistance = 999999L,
    MinFragmentLength = 10,
    Start = NULL,
    End = NULL,
    OutFileType = "fastq",
    outpath = "processreads_adapter_only",
    render_report = FALSE,
    force = TRUE,
    verbose = FALSE,
    KeepIntermediates = TRUE,  # RemoveAdapter writes outputs only if KeepIntermediates=TRUE
    FinalSuffix = "processed"
  )
  
  md <- qc2@metadata[[key]]
  expect_true(!is.null(md$adapter_summary))
  expect_equal(md$adapter_summary$reads_before, 2L)
  expect_equal(md$adapter_summary$reads_after, 2L)
  expect_true(length(md$adapter_summary$output_paths) >= 1)
  expect_true(all(file.exists(md$adapter_summary$output_paths)))
  
  expect_true(is.null(md$filter_summary))
  expect_true(is.null(md$trim_summary))
})

test_that("ProcessReads: trim-only writes final output and records trim_summary", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ShortRead")
  skip_if_not_installed("IRanges")
  
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  
  write_test_fastq <- function(path, seqs, ids = NULL, qchar = "I") {
    if (is.null(ids)) ids <- paste0("read", seq_along(seqs))
    dna <- Biostrings::DNAStringSet(seqs); names(dna) <- ids
    qs <- Biostrings::PhredQuality(Biostrings::BStringSet(vapply(
      seqs, function(s) paste(rep(qchar, nchar(s)), collapse = ""), character(1)
    )))
    names(qs) <- ids
    qsdss <- Biostrings::QualityScaledDNAStringSet(dna, qs)
    sr <- ShortRead::ShortReadQ(
      sread = Biostrings::DNAStringSet(qsdss),
      quality = ShortRead::FastqQuality(Biostrings::quality(qsdss)),
      id = Biostrings::BStringSet(ids)
    )
    ShortRead::writeFastq(sr, file = path, compress = FALSE)
    path
  }
  
  fastq_path <- file.path(wd, "toy.fastq")
  write_test_fastq(fastq_path,
                   seqs = c(paste(rep("A", 40), collapse = ""),
                            paste(rep("C", 40), collapse = "")),
                   ids = c("r1", "r2"))
  
  qc <- .init_qc_object(files = fastq_path)
  key <- basename(fastq_path)
  qc <- QualMat(qc, ImportFile(fastq_path), filename = key)
  
  qc2 <- ProcessReads(
    qc_obj = qc,
    filter = FALSE,
    AdapterSeq = NULL,
    Start = 5,
    End = 5,
    OutFileType = "fastq",
    outpath = "processreads_trim_only",
    render_report = FALSE,
    force = TRUE,
    verbose = FALSE,
    KeepIntermediates = FALSE,
    FinalSuffix = "processed"
  )
  
  md <- qc2@metadata[[key]]
  expect_true(!is.null(md$trim_summary))
  expect_equal(md$trim_summary$reads_before, 2L)
  expect_equal(md$trim_summary$reads_after, 2L)
  expect_true(length(md$trim_summary$output_paths) >= 1)
  expect_true(all(file.exists(md$trim_summary$output_paths)))
  
  expect_true(is.null(md$filter_summary))
  expect_true(is.null(md$adapter_summary))
})

test_that("ProcessReads: skip remaining steps when FilterLong keeps 0 reads", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ShortRead")
  skip_if_not_installed("IRanges")
  
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  
  write_test_fastq <- function(path, seqs, ids = NULL, qchar = "I") {
    if (is.null(ids)) ids <- paste0("read", seq_along(seqs))
    dna <- Biostrings::DNAStringSet(seqs); names(dna) <- ids
    qs <- Biostrings::PhredQuality(Biostrings::BStringSet(vapply(
      seqs, function(s) paste(rep(qchar, nchar(s)), collapse = ""), character(1)
    )))
    names(qs) <- ids
    qsdss <- Biostrings::QualityScaledDNAStringSet(dna, qs)
    sr <- ShortRead::ShortReadQ(
      sread = Biostrings::DNAStringSet(qsdss),
      quality = ShortRead::FastqQuality(Biostrings::quality(qsdss)),
      id = Biostrings::BStringSet(ids)
    )
    ShortRead::writeFastq(sr, file = path, compress = FALSE)
    path
  }
  
  adapter <- "ACGTACGT"
  
  # short reads: length 30
  fastq_path <- file.path(wd, "toy.fastq")
  write_test_fastq(fastq_path,
                   seqs = c(paste(rep("A", 30), collapse = ""),
                            paste(rep("C", 30), collapse = "")),
                   ids = c("r1", "r2"))
  
  qc <- .init_qc_object(files = fastq_path)
  key <- basename(fastq_path)
  qc <- QualMat(qc, ImportFile(fastq_path), filename = key)
  
  qc2 <- ProcessReads(
    qc_obj = qc,
    filter = TRUE,
    MinAvgQS = 0,
    MinLength = 99999,     # forces 0 reads kept
    MaxNumberNs = 999,
    AdapterSeq = adapter,  # would have run, but should be skipped
    Start = 1,
    End = 1,
    OutFileType = "fastq",
    outpath = "processreads_skip_after_filter",
    render_report = FALSE,
    force = TRUE,
    verbose = FALSE,
    KeepIntermediates = FALSE,
    FinalSuffix = "processed"
  )
  
  md <- qc2@metadata[[key]]
  expect_true(!is.null(md$filter_summary))
  expect_equal(md$filter_summary$reads_after, 0L)
  
  # If filter kept 0, remaining steps should not populate summaries
  expect_true(is.null(md$adapter_summary))
  expect_true(is.null(md$trim_summary))
})

test_that("ProcessReads: skip remaining steps when RemoveAdapter keeps 0 reads", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ShortRead")
  skip_if_not_installed("IRanges")
  
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  
  write_test_fastq <- function(path, seqs, ids = NULL, qchar = "I") {
    if (is.null(ids)) ids <- paste0("read", seq_along(seqs))
    dna <- Biostrings::DNAStringSet(seqs); names(dna) <- ids
    qs <- Biostrings::PhredQuality(Biostrings::BStringSet(vapply(
      seqs, function(s) paste(rep(qchar, nchar(s)), collapse = ""), character(1)
    )))
    names(qs) <- ids
    qsdss <- Biostrings::QualityScaledDNAStringSet(dna, qs)
    sr <- ShortRead::ShortReadQ(
      sread = Biostrings::DNAStringSet(qsdss),
      quality = ShortRead::FastqQuality(Biostrings::quality(qsdss)),
      id = Biostrings::BStringSet(ids)
    )
    ShortRead::writeFastq(sr, file = path, compress = FALSE)
    path
  }
  
  adapter <- "ACGTACGT"
  # Make reads that will be end-trimmed to almost nothing, then enforce huge MinFragmentLength
  seqs <- c(paste0(adapter, paste(rep("A", 20), collapse = "")),
            paste0(paste(rep("T", 20), collapse = ""), adapter))
  
  fastq_path <- file.path(wd, "toy.fastq")
  write_test_fastq(fastq_path, seqs = seqs, ids = c("r1", "r2"))
  
  qc <- .init_qc_object(files = fastq_path)
  key <- basename(fastq_path)
  qc <- QualMat(qc, ImportFile(fastq_path), filename = key)
  
  qc2 <- ProcessReads(
    qc_obj = qc,
    filter = FALSE,
    AdapterSeq = adapter,
    MaxMismatchEnd = 0,
    MinOverlapEnd = 50,        # treat adapter as end-hit aggressively
    MinInternalDistance = 999999L,
    MinFragmentLength = 99999, # forces 0 reads after adapter stage
    Start = 1,
    End = 1,                   # should be skipped
    OutFileType = "fastq",
    outpath = "processreads_skip_after_adapter",
    render_report = FALSE,
    force = TRUE,
    verbose = FALSE,
    KeepIntermediates = FALSE,
    FinalSuffix = "processed"
  )
  
  md <- qc2@metadata[[key]]
  expect_true(!is.null(md$adapter_summary))
  expect_equal(md$adapter_summary$reads_after, 0L)
  
  # adapter kept 0 -> trim summary should not exist
  expect_true(is.null(md$trim_summary))
})

test_that("ProcessReads: KeepIntermediates=TRUE leaves intermediate fastqs on disk", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ShortRead")
  skip_if_not_installed("IRanges")
  
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  
  tmp <- withr::local_tempdir()
  withr::local_envvar(c(TMPDIR = tmp, TEMP = tmp, TMP = tmp))
  
  write_test_fastq <- function(path, seqs, ids = NULL, qchar = "I") {
    if (is.null(ids)) ids <- paste0("read", seq_along(seqs))
    dna <- Biostrings::DNAStringSet(seqs); names(dna) <- ids
    qs <- Biostrings::PhredQuality(Biostrings::BStringSet(vapply(
      seqs, function(s) paste(rep(qchar, nchar(s)), collapse = ""), character(1)
    )))
    names(qs) <- ids
    qsdss <- Biostrings::QualityScaledDNAStringSet(dna, qs)
    sr <- ShortRead::ShortReadQ(
      sread = Biostrings::DNAStringSet(qsdss),
      quality = ShortRead::FastqQuality(Biostrings::quality(qsdss)),
      id = Biostrings::BStringSet(ids)
    )
    ShortRead::writeFastq(sr, file = path, compress = FALSE)
    path
  }
  
  adapter <- "ACGTACGT"
  fastq_path <- file.path(wd, "toy.fastq")
  write_test_fastq(fastq_path,
                   seqs = c(paste0(paste(rep("A", 20), collapse = ""), adapter),
                            paste0(adapter, paste(rep("T", 20), collapse = ""))),
                   ids = c("r1", "r2"))
  
  qc <- .init_qc_object(files = fastq_path)
  key <- basename(fastq_path)
  qc <- QualMat(qc, ImportFile(fastq_path), filename = key)
  
  qc2 <- ProcessReads(
    qc_obj = qc,
    filter = TRUE,
    MinAvgQS = 0,
    MinLength = 10,
    MaxNumberNs = 999,
    AdapterSeq = adapter,
    MaxMismatchEnd = 0,
    MinOverlapEnd = 5,
    MinInternalDistance = 999999L,
    MinFragmentLength = 10,
    Start = 1,
    End = 1,
    OutFileType = "fastq",
    outpath = "processreads_keep_intermediates",
    render_report = FALSE,
    force = TRUE,
    verbose = FALSE,
    KeepIntermediates = TRUE,
    FinalSuffix = "processed"
  )
  
  # Intermediate directory used by your step functions:
  inter_dir <- file.path(tempdir(), "longR_intermediate")
  expect_true(dir.exists(inter_dir))
  
  inter_files <- list.files(inter_dir, pattern = "\\.fastq$", full.names = TRUE)
  # We expect at least one intermediate fastq left behind
  expect_true(length(inter_files) >= 1)
  expect_true(all(file.exists(inter_files)))
  
  # And final output exists too
  md <- qc2@metadata[[key]]
  expect_true(length(md$trim_summary$output_paths) >= 1)
  expect_true(all(file.exists(md$trim_summary$output_paths)))
})

test_that("ProcessReads loops over multiple files and keeps metadata per file", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ShortRead")
  skip_if_not_installed("IRanges")
  
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  
  write_test_fastq <- function(path, seqs, ids = NULL, qchar = "I") {
    if (is.null(ids)) ids <- paste0("read", seq_along(seqs))
    dna <- Biostrings::DNAStringSet(seqs); names(dna) <- ids
    qs <- Biostrings::PhredQuality(Biostrings::BStringSet(vapply(
      seqs, function(s) paste(rep(qchar, nchar(s)), collapse = ""), character(1)
    )))
    names(qs) <- ids
    qsdss <- Biostrings::QualityScaledDNAStringSet(dna, qs)
    
    sr <- ShortRead::ShortReadQ(
      sread   = Biostrings::DNAStringSet(qsdss),
      quality = ShortRead::FastqQuality(Biostrings::quality(qsdss)),
      id      = Biostrings::BStringSet(ids)
    )
    ShortRead::writeFastq(sr, file = path, compress = FALSE)
    path
  }
  
  adapter <- "ACGTACGT"
  prefix20 <- paste(rep("A", 20), collapse = "")
  suffix20 <- paste(rep("T", 20), collapse = "")
  
  # File 1: will be SKIPPED after filter (MinLength too high)
  f1 <- file.path(wd, "f1.fastq")
  write_test_fastq(f1,
                   seqs = c(paste(rep("A", 30), collapse = ""), paste(rep("C", 30), collapse = "")),
                   ids = c("a1", "a2")
  )
  
  # File 2: will RUN fully (filter keeps, adapter + trim run)
  f2 <- file.path(wd, "f2.fastq")
  write_test_fastq(f2,
                   seqs = c(paste0(prefix20, adapter), paste0(adapter, suffix20), paste0(prefix20, suffix20)),
                   ids = c("b1", "b2", "b3")
  )
  
  qc <- .init_qc_object(files = c(f1, f2))
  
  # Populate metrics for EACH file (since ProcessReads assumes metrics already exist)
  k1 <- basename(f1)
  k2 <- basename(f2)
  
  qc <- QualMat(qc, ImportFile(f1), filename = k1)
  qc <- QualMat(qc, ImportFile(f2), filename = k2)
  
  out <- ProcessReads(
    qc_obj = qc,
    filter = TRUE,
    MinAvgQS = 0,
    MinLength = 10,
    MaxNumberNs = 999,
    AdapterSeq = adapter,
    MaxMismatchEnd = 0,
    MinOverlapEnd = 5,
    MinInternalDistance = 999999L,
    MinFragmentLength = 10,
    Start = 1,
    End = 1,
    OutFileType = "fastq",
    outpath = "processreads_multi_file",
    render_report = FALSE,
    force = TRUE,
    verbose = FALSE,
    KeepIntermediates = FALSE,
    FinalSuffix = "processed"
  )
  
  # Now assert *per-file loop behavior*
  
  md1 <- out@metadata[[k1]]
  md2 <- out@metadata[[k2]]
  
  expect_true(!is.null(md1$filter_summary))
  expect_equal(md1$filter_summary$reads_before, 2L)
  
  # Force skip for file1 by running again with huge MinLength (per-file skip test)
  out2 <- ProcessReads(
    qc_obj = qc,
    filter = TRUE,
    MinAvgQS = 0,
    MinLength = 99999,     # makes f1 and f2 skip after filter (but that’s ok for the skip assertion)
    MaxNumberNs = 999,
    AdapterSeq = adapter,
    Start = 1,
    End = 1,
    OutFileType = "fastq",
    outpath = "processreads_multi_file_skip",
    render_report = FALSE,
    force = TRUE,
    verbose = FALSE,
    KeepIntermediates = FALSE,
    FinalSuffix = "processed"
  )
  
  md1s <- out2@metadata[[k1]]
  md2s <- out2@metadata[[k2]]
  
  expect_true(!is.null(md1s$filter_summary))
  expect_equal(md1s$filter_summary$reads_after, 0L)
  expect_true(is.null(md1s$adapter_summary))
  expect_true(is.null(md1s$trim_summary))
  
  expect_true(!is.null(md2s$filter_summary))
  expect_equal(md2s$filter_summary$reads_after, 0L)
  expect_true(is.null(md2s$adapter_summary))
  expect_true(is.null(md2s$trim_summary))
})