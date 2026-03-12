library(testthat)
library(WeinR)

test_that("CreateReport errors when directory has no supported files", {
  empty_dir <- tempfile("empty_dir_")
  dir.create(empty_dir)
  
  expect_error(
    CreateReport(
      files = empty_dir,
      report_name = "test_report",
      render_report = FALSE,
      force = TRUE
    ),
    "No .fastq, .fastq.gz, or .bam files found"
  )
})

test_that("CreateReport expands a directory of supported files", {
  input_dir <- tempfile("input_dir_")
  dir.create(input_dir)
  
  file.create(file.path(input_dir, "a.fastq"))
  file.create(file.path(input_dir, "b.fastq.gz"))
  file.create(file.path(input_dir, "c.bam"))
  file.create(file.path(input_dir, "ignore.txt"))
  
  captured_files <- NULL
  
  fake_qc <- methods::new(
    "LongReadQC",
    files = character(),
    summary_metrics = data.frame(File = character(), stringsAsFactors = FALSE),
    metrics = list()
  )
  
  local_mocked_bindings(
    .init_qc_object = function(files) {
      captured_files <<- sort(basename(files))
      fake_qc@files <- files
      fake_qc
    },
    ImportFile = function(file, MinNumReads) {
      list(file = file, MinNumReads = MinNumReads)
    },
    QualMat = function(qc_obj, qsds, fname) {
      qc_obj@metrics[[fname]] <- list(
        metric1 = data.frame(x = 1)
      )
      qc_obj@summary_metrics <- rbind(
        qc_obj@summary_metrics,
        data.frame(File = fname, stringsAsFactors = FALSE)
      )
      qc_obj
    },
    QualPlot = function(qc_obj, filename, out_dir, dpi, overwrite) {
      qc_obj
    },
    Reporter = function(...) {
      invisible(list(rmd = "fake.Rmd", html = "fake.html"))
    }
  )
  
  old_wd <- getwd()
  tmp_wd <- tempfile("wd_")
  dir.create(tmp_wd)
  setwd(tmp_wd)
  withr::defer(setwd(old_wd))
  
  res <- CreateReport(
    files = input_dir,
    report_name = "dir_report",
    render_report = FALSE,
    force = TRUE
  )
  
  expect_s4_class(res, "LongReadQC")
  expect_equal(captured_files, sort(c("a.fastq", "b.fastq.gz", "c.bam")))
})

test_that("CreateReport returns a LongReadQC object and writes outputs", {
  fake_files <- c("sample1.fastq", "sample2.fastq")
  
  fake_qc <- methods::new(
    "LongReadQC",
    files = fake_files,
    summary_metrics = data.frame(File = character(), stringsAsFactors = FALSE),
    metrics = list()
  )
  
  local_mocked_bindings(
    .init_qc_object = function(files) {
      fake_qc@files <- files
      fake_qc
    },
    ImportFile = function(file, MinNumReads) {
      list(file = file, MinNumReads = MinNumReads)
    },
    QualMat = function(qc_obj, qsds, fname) {
      qc_obj@metrics[[fname]] <- list(
        read_lengths = data.frame(length = c(100, 200, 300))
      )
      qc_obj@summary_metrics <- rbind(
        qc_obj@summary_metrics,
        data.frame(File = fname, stringsAsFactors = FALSE)
      )
      qc_obj
    },
    QualPlot = function(qc_obj, filename, out_dir, dpi, overwrite) {
      qc_obj
    },
    Reporter = function(...) {
      invisible(list(rmd = "fake.Rmd", html = "fake.html"))
    }
  )
  
  old_wd <- getwd()
  tmp_wd <- tempfile("wd_")
  dir.create(tmp_wd)
  setwd(tmp_wd)
  withr::defer(setwd(old_wd))
  
  res <- CreateReport(
    files = fake_files,
    report_name = "my_report",
    render_report = FALSE,
    force = TRUE
  )
  
  expect_s4_class(res, "LongReadQC")
  
  out_root <- file.path(tmp_wd, "WeinR_Outputs")
  expect_true(dir.exists(out_root))
  
  run_dirs <- list.dirs(out_root, recursive = FALSE, full.names = TRUE)
  expect_length(run_dirs, 1)
  
  run_dir <- run_dirs[[1]]
  
  expect_true(dir.exists(file.path(run_dir, "metrics")))
  expect_true(dir.exists(file.path(run_dir, "plots")))
  expect_true(dir.exists(file.path(run_dir, "reports")))
  expect_true(dir.exists(file.path(run_dir, "objects")))
  
  expect_true(file.exists(file.path(run_dir, "objects", "qc.rds")))
  expect_true(file.exists(file.path(run_dir, "metrics", "my_report_summary_metrics.csv")))
  expect_true(file.exists(file.path(run_dir, "metrics", "sample1.fastq", "read_lengths.csv")))
  expect_true(file.exists(file.path(run_dir, "metrics", "sample2.fastq", "read_lengths.csv")))
})

test_that("CreateReport calls Reporter only when render_report is TRUE", {
  fake_files <- "sample1.fastq"
  reporter_called <- FALSE
  
  fake_qc <- methods::new(
    "LongReadQC",
    files = fake_files,
    summary_metrics = data.frame(File = character(), stringsAsFactors = FALSE),
    metrics = list()
  )
  
  local_mocked_bindings(
    .init_qc_object = function(files) {
      fake_qc@files <- files
      fake_qc
    },
    ImportFile = function(file, MinNumReads) {
      list(file = file)
    },
    QualMat = function(qc_obj, qsds, fname) {
      qc_obj@metrics[[fname]] <- list(
        metric1 = data.frame(x = 1)
      )
      qc_obj@summary_metrics <- rbind(
        qc_obj@summary_metrics,
        data.frame(File = fname, stringsAsFactors = FALSE)
      )
      qc_obj
    },
    QualPlot = function(qc_obj, filename, out_dir, dpi, overwrite) {
      qc_obj
    },
    Reporter = function(...) {
      reporter_called <<- TRUE
      invisible(list(rmd = "fake.Rmd", html = "fake.html"))
    }
  )
  
  old_wd <- getwd()
  tmp_wd <- tempfile("wd_")
  dir.create(tmp_wd)
  setwd(tmp_wd)
  withr::defer(setwd(old_wd))
  
  CreateReport(
    files = fake_files,
    report_name = "report_true",
    render_report = TRUE,
    force = TRUE
  )
  
  expect_true(reporter_called)
  
  reporter_called <- FALSE
  
  CreateReport(
    files = fake_files,
    report_name = "report_false",
    render_report = FALSE,
    force = TRUE
  )
  
  expect_false(reporter_called)
})

test_that("CreateReport errors if summary_metrics lacks File column", {
  fake_files <- "sample1.fastq"
  
  fake_qc <- methods::new(
    "LongReadQC",
    files = fake_files,
    summary_metrics = data.frame(NotFile = 1),
    metrics = list()
  )
  
  local_mocked_bindings(
    .init_qc_object = function(files) {
      fake_qc@files <- files
      fake_qc
    },
    ImportFile = function(file, MinNumReads) {
      list(file = file)
    },
    QualMat = function(qc_obj, qsds, fname) {
      qc_obj@metrics[[fname]] <- list(
        metric1 = data.frame(x = 1)
      )
      qc_obj
    },
    QualPlot = function(qc_obj, filename, out_dir, dpi, overwrite) {
      qc_obj
    },
    Reporter = function(...) {
      invisible(list(rmd = "fake.Rmd", html = "fake.html"))
    }
  )
  
  old_wd <- getwd()
  tmp_wd <- tempfile("wd_")
  dir.create(tmp_wd)
  setwd(tmp_wd)
  withr::defer(setwd(old_wd))
  
  expect_error(
    CreateReport(
      files = fake_files,
      report_name = "bad_summary",
      render_report = FALSE,
      force = TRUE
    ),
    "must contain a 'File' column"
  )
})

test_that("CreateReport skips files when ImportFile returns NULL", {
  fake_files <- c("keep.fastq", "skip.fastq")
  
  fake_qc <- methods::new(
    "LongReadQC",
    files = fake_files,
    summary_metrics = data.frame(File = character(), stringsAsFactors = FALSE),
    metrics = list()
  )
  
  local_mocked_bindings(
    .init_qc_object = function(files) {
      fake_qc@files <- files
      fake_qc
    },
    ImportFile = function(file, MinNumReads) {
      if (basename(file) == "skip.fastq") return(NULL)
      list(file = file)
    },
    QualMat = function(qc_obj, qsds, fname) {
      qc_obj@metrics[[fname]] <- list(
        metric1 = data.frame(x = 1)
      )
      qc_obj@summary_metrics <- rbind(
        qc_obj@summary_metrics,
        data.frame(File = fname, stringsAsFactors = FALSE)
      )
      qc_obj
    },
    QualPlot = function(qc_obj, filename, out_dir, dpi, overwrite) {
      qc_obj
    },
    Reporter = function(...) {
      invisible(list(rmd = "fake.Rmd", html = "fake.html"))
    }
  )
  
  old_wd <- getwd()
  tmp_wd <- tempfile("wd_")
  dir.create(tmp_wd)
  setwd(tmp_wd)
  withr::defer(setwd(old_wd))
  
  res <- CreateReport(
    files = fake_files,
    report_name = "skip_test",
    render_report = FALSE,
    force = TRUE
  )
  
  expect_true("keep.fastq" %in% names(res@metrics))
  expect_false("skip.fastq" %in% names(res@metrics))
  expect_equal(res@summary_metrics$File, "keep.fastq")
})