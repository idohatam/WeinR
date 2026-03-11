# tests/testthat/test-FastqCheck.R

test_that("FastqCheck validates n_records", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "x.fastq")
  file.create(f)
  
  expect_error(FastqCheck(f, n_records = NULL), "n_records must be a single positive integer")
  expect_error(FastqCheck(f, n_records = "10"), "n_records must be a single positive integer")
  expect_error(FastqCheck(f, n_records = c(1,2)), "n_records must be a single positive integer")
  expect_error(FastqCheck(f, n_records = NA), "n_records must be a single positive integer")
  expect_error(FastqCheck(f, n_records = 0), "n_records must be a single positive integer")
  expect_error(FastqCheck(f, n_records = 1.5), "n_records must be a single positive integer")
})

test_that("FastqCheck treats empty file as valid", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "empty.fastq")
  file.create(f)
  expect_true(isTRUE(FastqCheck(f, n_records = 5)))
})

test_that("FastqCheck errors on incomplete records (not multiple of 4)", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "bad.fastq")
  writeLines(c("@r1","ACGT","+"), f) # 3 lines
  expect_error(FastqCheck(f, n_records = 10),
               "FASTQ appears malformed: incomplete record \\(lines not multiple of 4\\)\\.")
})

test_that("FastqCheck errors if header lines don't start with @", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "bad_header.fastq")
  writeLines(c("r1","ACGT","+","FFFF"), f)
  expect_error(FastqCheck(f, n_records = 10),
               "FASTQ appears malformed: header line\\(s\\) not starting with '@'\\.")
})

test_that("FastqCheck errors if plus lines don't start with +", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "bad_plus.fastq")
  writeLines(c("@r1","ACGT","plus","FFFF"), f)
  expect_error(FastqCheck(f, n_records = 10),
               "FASTQ appears malformed: plus line\\(s\\) not starting with '\\+'\\.")
})

test_that("FastqCheck works on .gz inputs too", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "ok.fastq.gz")
  
  write_fastq_gz(f, minimal_records_ok(2))
  
  expect_true(isTRUE(FastqCheck(f, n_records = 10)))
})

