test_that("CheckFile validates inputs (type/length/NA/empty)", {
  expect_error(CheckFile(NULL), "File path is NULL\\.")
  expect_error(CheckFile(123), "File path must be a character string\\.")
  expect_error(CheckFile(c("a","b")), "File path must be a single path \\(length 1\\)\\.")
  expect_error(CheckFile(NA_character_), "File path is NA\\.")
  expect_error(CheckFile("   "), "File path is empty\\.")
})

test_that("CheckFile errors for nonexistent paths and directories", {
  p <- tempfile()
  expect_error(CheckFile(p), "File does not exist:")
  
  d <- tempfile()
  dir.create(d)
  expect_error(CheckFile(d), "Expected a file path, got a directory:")
})

test_that("CheckFile infers file types correctly (case-insensitive)", {
  tmp <- withr::local_tempdir()
  
  bam <- file.path(tmp, "x.BAM")
  fq  <- file.path(tmp, "x.Fq")
  fqgz <- file.path(tmp, "x.fastq.GZ")
  
  file.create(bam)
  file.create(fq)
  file.create(fqgz)
  
  expect_identical(CheckFile(bam), "bam")
  expect_identical(CheckFile(fq), "fastq")
  expect_identical(CheckFile(fqgz), "fastq.gz")
})

test_that("CheckFile errors on unsupported extensions", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "x.txt")
  file.create(f)
  expect_error(CheckFile(f), "Unsupported file type\\. Expected \\.bam or \\.fastq/\\.fq \\(optionally \\.gz\\)\\.")
})

## Didn't test the warning for larger files >=20 gb