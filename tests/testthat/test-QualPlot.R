library(testthat)
library(WeinR)

# ---- Helpers ----

# Builds a minimal valid LongReadQC object with one file entry
make_fake_qc <- function(filename = "sample.fastq") {
  obj <- new("LongReadQC",
             files = filename,
             metrics = list(),
             plots = list(),
             summary_metrics = data.frame(),
             metadata = list()
  )
  
  obj@metrics[[filename]] <- list(
    readLengths   = as.list(c(500L, 1000L, 1500L, 2000L, 3000L)),
    meanprQscore  = as.list(c(12.5, 15.0, 18.3, 20.1, 22.4)),
    prGCcontent   = as.list(c(0.45, 0.50, 0.55, 0.48, 0.52)),
    perPosQuality = data.frame(
      position = 1:10,
      mean     = rep(18, 10),
      median   = rep(17, 10),
      q25      = rep(14, 10),
      q75      = rep(21, 10)
    )
  )
  
  obj
}

# ---- Input validation ----

test_that("QualPlot errors when filename is not a character", {
  obj <- make_fake_qc()
  expect_error(QualPlot(obj, filename = 123), "`filename` must be a length-1 character string")
})

test_that("QualPlot errors when filename has length > 1", {
  obj <- make_fake_qc()
  expect_error(
    QualPlot(obj, filename = c("a.fastq", "b.fastq")),
    "`filename` must be a length-1 character string"
  )
})

test_that("QualPlot errors when filename is not in qc_obj@metrics", {
  obj <- make_fake_qc("sample.fastq")
  expect_error(
    QualPlot(obj, filename = "nonexistent.fastq"),
    "not found in qc_obj@metrics"
  )
})

test_that("QualPlot errors when metrics slot is empty and filename missing", {
  obj <- new("LongReadQC",
             files = character(0),
             metrics = list(),
             plots = list(),
             summary_metrics = data.frame(),
             metadata = list()
  )
  expect_error(
    QualPlot(obj, filename = "sample.fastq"),
    "not found in qc_obj@metrics"
  )
})

test_that("QualPlot errors when a required metric is missing", {
  obj <- make_fake_qc("sample.fastq")
  # Remove a required metric
  obj@metrics[["sample.fastq"]]$perPosQuality <- NULL
  expect_error(
    QualPlot(obj, filename = "sample.fastq"),
    "Missing metrics"
  )
})

test_that("QualPlot errors when only readLengths is missing", {
  obj <- make_fake_qc("sample.fastq")
  obj@metrics[["sample.fastq"]]$readLengths <- NULL
  expect_error(
    QualPlot(obj, filename = "sample.fastq"),
    "readLengths"
  )
})

test_that("QualPlot errors when only meanprQscore is missing", {
  obj <- make_fake_qc("sample.fastq")
  obj@metrics[["sample.fastq"]]$meanprQscore <- NULL
  expect_error(
    QualPlot(obj, filename = "sample.fastq"),
    "meanprQscore"
  )
})

test_that("QualPlot errors when only prGCcontent is missing", {
  obj <- make_fake_qc("sample.fastq")
  obj@metrics[["sample.fastq"]]$prGCcontent <- NULL
  expect_error(
    QualPlot(obj, filename = "sample.fastq"),
    "prGCcontent"
  )
})

test_that("QualPlot errors when only perPosQuality is missing", {
  obj <- make_fake_qc("sample.fastq")
  obj@metrics[["sample.fastq"]]$perPosQuality <- NULL
  expect_error(
    QualPlot(obj, filename = "sample.fastq"),
    "perPosQuality"
  )
})

test_that("QualPlot error message lists all missing metrics when multiple are absent", {
  obj <- make_fake_qc("sample.fastq")
  obj@metrics[["sample.fastq"]]$readLengths  <- NULL
  obj@metrics[["sample.fastq"]]$meanprQscore <- NULL
  err <- expect_error(QualPlot(obj, filename = "sample.fastq"))
  expect_match(err$message, "readLengths")
  expect_match(err$message, "meanprQscore")
})

# ---- Normal operation ----

test_that("QualPlot returns a LongReadQC object", {
  obj <- make_fake_qc("sample.fastq")
  out_dir <- tempfile("plots")
  result <- QualPlot(obj, filename = "sample.fastq", out_dir = out_dir)
  expect_s4_class(result, "LongReadQC")
})

test_that("QualPlot populates plots slot for the given filename", {
  obj <- make_fake_qc("sample.fastq")
  out_dir <- tempfile("plots")
  result <- QualPlot(obj, filename = "sample.fastq", out_dir = out_dir)
  expect_true("sample.fastq" %in% names(result@plots))
})

test_that("QualPlot stores all five expected plot keys", {
  obj <- make_fake_qc("sample.fastq")
  out_dir <- tempfile("plots")
  result <- QualPlot(obj, filename = "sample.fastq", out_dir = out_dir)
  plot_keys <- names(result@plots[["sample.fastq"]])
  expect_setequal(plot_keys, c("length_hist", "quality_hist", "gc_hist", "q_vs_length", "per_pos_q"))
})

test_that("QualPlot writes PNG files to disk", {
  obj <- make_fake_qc("sample.fastq")
  out_dir <- tempfile("plots")
  result <- QualPlot(obj, filename = "sample.fastq", out_dir = out_dir)
  paths <- unlist(result@plots[["sample.fastq"]])
  real_paths <- paths[!is.na(paths)]
  expect_true(all(file.exists(real_paths)))
})

test_that("QualPlot creates output directory if it does not exist", {
  obj <- make_fake_qc("sample.fastq")
  out_dir <- file.path(tempdir(), paste0("brand_new_dir_", as.integer(Sys.time()), sample(1e6, 1)))
  expect_false(dir.exists(out_dir))
  QualPlot(obj, filename = "sample.fastq", out_dir = out_dir)
  expect_true(dir.exists(out_dir))
})

test_that("QualPlot does not overwrite existing plots when overwrite = FALSE", {
  obj <- make_fake_qc("sample.fastq")
  out_dir <- tempfile("plots")
  
  # First run — creates the files
  result1 <- QualPlot(obj, filename = "sample.fastq", out_dir = out_dir, overwrite = FALSE)
  paths <- unlist(result1@plots[["sample.fastq"]])
  real_paths <- paths[!is.na(paths)]
  mtimes_before <- file.mtime(real_paths)
  
  Sys.sleep(1)
  
  # Second run — should not touch existing files
  QualPlot(obj, filename = "sample.fastq", out_dir = out_dir, overwrite = FALSE)
  mtimes_after <- file.mtime(real_paths)
  
  expect_equal(mtimes_before, mtimes_after)
})

test_that("QualPlot sanitizes filenames with special characters in directory name", {
  fancy_name <- "path/to/my file (1).fastq"
  obj <- make_fake_qc(fancy_name)
  out_dir <- tempfile("plots")
  # Should not error on special characters
  expect_no_error(QualPlot(obj, filename = fancy_name, out_dir = out_dir))
})

test_that("QualPlot preserves pre-existing plots for other files", {
  obj <- make_fake_qc("sample.fastq")
  obj@plots[["other.fastq"]] <- list(length_hist = "fake/path.png")
  out_dir <- tempfile("plots")
  result <- QualPlot(obj, filename = "sample.fastq", out_dir = out_dir)
  expect_true("other.fastq" %in% names(result@plots))
})

test_that("QualPlot stores NA for plots that return NULL (empty GC content)", {
  obj <- make_fake_qc("sample.fastq")
  obj@metrics[["sample.fastq"]]$prGCcontent <- list()  # empty -> make_gc_plot returns NULL
  out_dir <- tempfile("plots")
  result <- QualPlot(obj, filename = "sample.fastq", out_dir = out_dir)
  expect_true(is.na(result@plots[["sample.fastq"]]$gc_hist))
})

test_that("QualPlot errors when a file's metrics entry is not a list", {
  obj <- make_fake_qc("sample.fastq")
  obj@metrics[["sample.fastq"]] <- "not a list"
  expect_error(
    QualPlot(obj, filename = "sample.fastq"),
    "must be a non-NULL list"
  )
})

# ---- Integration test with real ONT data ----

test_that("QualPlot works end-to-end on real ONT data", {
  file <- system.file("extdata", "sample.fastq", package = "WeinR")
  skip_if(file == "", "sample.fastq not available")
  
  qc_obj <- WeinR:::.init_qc_object(file)
  stringset <- WeinR:::ImportFile(file)
  qc_obj <- WeinR:::QualMat(qc_obj, stringset, file)
  
  out_dir <- tempfile("plots")
  result <- QualPlot(qc_obj, filename = file, out_dir = out_dir)
  
  expect_s4_class(result, "LongReadQC")
  expect_true(file %in% names(result@plots))
  expect_length(result@plots[[file]], 5)
  paths <- unlist(result@plots[[file]])
  real_paths <- paths[!is.na(paths)]
  expect_true(all(file.exists(real_paths)))
})