test_that("ImportFile validates MinNumReads", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "ok.fastq")
  write_fastq(f, minimal_records_ok(1))
  
  expect_error(ImportFile(f, MinNumReads = NULL), "MinNumReads must be a single positive integer")
  expect_error(ImportFile(f, MinNumReads = "2"), "MinNumReads must be a single positive integer")
  expect_error(ImportFile(f, MinNumReads = c(1,2)), "MinNumReads must be a single positive integer")
  expect_error(ImportFile(f, MinNumReads = NA), "MinNumReads must be a single positive integer")
  expect_error(ImportFile(f, MinNumReads = 0), "MinNumReads must be a single positive integer")
  expect_error(ImportFile(f, MinNumReads = 1.2), "MinNumReads must be a single positive integer")
})

test_that("ImportFile errors are propagated from CheckFile for bad paths/types", {
  expect_error(ImportFile(NULL), "File path is NULL\\.")
  expect_error(ImportFile("nope.fastq"), "File does not exist:")
})

test_that("ImportFile FASTQ branch: FastqCheck failure is rethrown with prefix", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "bad.fastq")
  
  # Structurally broken: 3 lines only
  writeLines(c("@r1","ACGT","+"), f)
  
  expect_error(
    ImportFile(f),
    "FASTQ check failed \\(corrupted/malformed\\): FASTQ appears malformed: incomplete record"
  )
})

test_that("ImportFile FASTQ branch: successful import returns QualityScaledDNAStringSet", {
  testthat::skip_if_not_installed("Biostrings")
  
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "ok.fastq")
  write_fastq(f, minimal_records_ok(3))
  
  out <- ImportFile(f)
  expect_s4_class(out, "QualityScaledDNAStringSet")
  expect_equal(length(out), 3L)
})

test_that("ImportFile FASTQ.GZ branch works", {
  testthat::skip_if_not_installed("Biostrings")
  
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "ok.fastq.gz")
  write_fastq_gz(f, minimal_records_ok(2))
  
  out <- ImportFile(f)
  expect_s4_class(out, "QualityScaledDNAStringSet")
  expect_equal(length(out), 2L)
})

test_that("ImportFile skips file with < MinNumReads and returns NULL with warning", {
  testthat::skip_if_not_installed("Biostrings")
  
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "ok.fastq")
  write_fastq(f, minimal_records_ok(2))
  
  expect_warning(
    out <- ImportFile(f, MinNumReads = 3),
    "File has 2 reads \\(< MinNumReads = 3\\)\\. File skipped:"
  )
  expect_null(out)
})


# bam files

test_that("ImportFile BAM branch works on Rsamtools extdata example", {
  testthat::skip_if_not_installed("Rsamtools")
  testthat::skip_if_not_installed("Biostrings")
  
  bam_path <- system.file("extdata", "ex1.bam", package = "Rsamtools")
  if (bam_path == "") testthat::skip("Rsamtools extdata BAM not found")
  
  out <- ImportFile(bam_path, MinNumReads = 1L)
  expect_s4_class(out, "QualityScaledDNAStringSet")
  expect_gt(length(out), 0L)
})

test_that("ImportFile BAM branch respects MinNumReads (warn + NULL)", {
  testthat::skip_if_not_installed("Rsamtools")
  testthat::skip_if_not_installed("Biostrings")
  
  bam_path <- system.file("extdata", "ex1.bam", package = "Rsamtools")
  if (bam_path == "") testthat::skip("Rsamtools extdata BAM not found")
  
  # Set threshold absurdly high so it must skip
  expect_warning(
    out <- ImportFile(bam_path, MinNumReads = 1e9),
    "File has .* reads \\(< MinNumReads = 1000000000\\)\\. File skipped:"
  )
  expect_null(out)
})

## Didn't check fastq files with different qual and seq lengths as Biostrings wasn't throwing an error for them. 