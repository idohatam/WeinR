library(testthat)
library(WeinR)

# ---- make_length_plot ----

test_that("make_length_plot returns a ggplot object", {
  lengths <- as.list(c(500, 1000, 2000, 5000, 10000))
  p <- make_length_plot(lengths, file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_length_plot accepts a plain numeric vector", {
  p <- make_length_plot(c(100, 200, 300), file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_length_plot silently drops zero and negative lengths", {
  lengths <- c(0, -5, 100, 500, 1000)
  # Should not error or warn — zeros filtered before log10
  expect_no_warning(make_length_plot(lengths, file_name = "test.fastq"))
})

test_that("make_length_plot silently drops non-finite values", {
  lengths <- c(NA, Inf, -Inf, NaN, 500, 1000)
  expect_no_error(make_length_plot(lengths, file_name = "test.fastq"))
})

test_that("make_length_plot handles a single valid read", {
  p <- make_length_plot(list(1000), file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_length_plot uses file_name as plot title", {
  p <- make_length_plot(c(500, 1000), file_name = "my_sample.fastq")
  expect_equal(p$labels$title, "my_sample.fastq")
})

# ---- make_quality_plot ----

test_that("make_quality_plot returns a ggplot object for valid input", {
  scores <- as.list(c(10, 15, 18, 20, 22, 25))
  p <- make_quality_plot(scores, file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_quality_plot returns NULL for empty input", {
  expect_null(make_quality_plot(list(), file_name = "test.fastq"))
})

test_that("make_quality_plot returns NULL when all values are non-finite", {
  expect_null(make_quality_plot(list(NA, NaN, Inf), file_name = "test.fastq"))
})

test_that("make_quality_plot handles a single read", {
  p <- make_quality_plot(list(20.0), file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_quality_plot uses provided avgQScore instead of computing mean", {
  scores <- as.list(c(10, 20, 30))
  p <- make_quality_plot(scores, avgQScore = 99, file_name = "test.fastq")
  # Extract the xintercept from geom_vline layers
  vline_data <- lapply(p$layers, function(l) l$geom_params)
  xintercepts <- sapply(p$layers, function(l) {
    if (inherits(l$geom, "GeomVline")) l$data$xintercept else NULL
  })
  xintercepts <- unlist(xintercepts)
  expect_true(99 %in% xintercepts)
})

test_that("make_quality_plot handles extreme Q-score values (0 and 40)", {
  scores <- as.list(c(0, 0, 40, 40))
  p <- make_quality_plot(scores, file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

# ---- make_gc_plot ----

test_that("make_gc_plot returns a ggplot object for valid input", {
  gc <- as.list(c(0.3, 0.45, 0.5, 0.55, 0.6))
  p <- make_gc_plot(gc, file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_gc_plot returns NULL for empty input", {
  expect_null(make_gc_plot(list(), file_name = "test.fastq"))
})

test_that("make_gc_plot returns NULL when all values are non-finite", {
  expect_null(make_gc_plot(list(NA, NaN, Inf), file_name = "test.fastq"))
})

test_that("make_gc_plot handles a single read", {
  p <- make_gc_plot(list(0.50), file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_gc_plot handles extreme GC values (0 and 1)", {
  # Reads that are entirely A/T or entirely G/C
  gc <- as.list(c(0.0, 0.0, 1.0, 1.0))
  p <- make_gc_plot(gc, file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_gc_plot uses file_name as plot title", {
  p <- make_gc_plot(list(0.5, 0.6), file_name = "my_sample.fastq")
  expect_equal(p$labels$title, "my_sample.fastq")
})

# ---- make_quality_vs_length_plot ----

test_that("make_quality_vs_length_plot returns a ggplot for valid input", {
  lengths <- as.list(c(500, 1000, 2000, 5000))
  scores  <- as.list(c(15,  18,   20,   22))
  p <- make_quality_vs_length_plot(lengths, scores, file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_quality_vs_length_plot returns NULL when all lengths are zero or negative", {
  lengths <- list(0, -1, -100)
  scores  <- list(15, 18, 20)
  expect_null(make_quality_vs_length_plot(lengths, scores, file_name = "test.fastq"))
})

test_that("make_quality_vs_length_plot returns NULL when all Q-scores are non-finite", {
  lengths <- list(500, 1000)
  scores  <- list(NA, NaN)
  expect_null(make_quality_vs_length_plot(lengths, scores, file_name = "test.fastq"))
})

test_that("make_quality_vs_length_plot handles a single valid point", {
  p <- make_quality_vs_length_plot(list(1000), list(20), file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_quality_vs_length_plot handles mixed valid and invalid points", {
  lengths <- list(0, NA, 500, 1000)
  scores  <- list(10, 20, 15,  18)
  # Should not error — invalid rows are filtered out
  expect_no_error(make_quality_vs_length_plot(lengths, scores, file_name = "test.fastq"))
})

# ---- make_per_position_plot ----

make_fake_per_pos <- function(n = 20) {
  data.frame(
    position = seq_len(n),
    mean     = rep(18, n),
    median   = rep(17, n),
    q25      = rep(14, n),
    q75      = rep(21, n)
  )
}

test_that("make_per_position_plot returns a ggplot for valid input", {
  p <- make_per_position_plot(make_fake_per_pos(), file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_per_position_plot returns NULL for NULL input", {
  expect_null(make_per_position_plot(NULL, file_name = "test.fastq"))
})

test_that("make_per_position_plot returns NULL for empty data frame", {
  empty_df <- data.frame(position = integer(0), mean = numeric(0),
                         median = numeric(0), q25 = numeric(0), q75 = numeric(0))
  expect_null(make_per_position_plot(empty_df, file_name = "test.fastq"))
})

test_that("make_per_position_plot errors when required columns are missing", {
  bad_df <- data.frame(position = 1:5, mean = rep(18, 5))  # missing median, q25, q75
  expect_error(
    make_per_position_plot(bad_df, file_name = "test.fastq"),
    "must contain columns"
  )
})

test_that("make_per_position_plot handles a single position row", {
  df <- data.frame(position = 1, mean = 18, median = 17, q25 = 14, q75 = 21)
  p <- make_per_position_plot(df, file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_per_position_plot skips rolling smoothing when nrow <= smoothing window", {
  # 10 rows is well below the 5001-row smoothing threshold — should still return a plot
  p <- make_per_position_plot(make_fake_per_pos(10), file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_per_position_plot applies rolling smoothing for large inputs", {
  # 6000 rows is above the 5001-row smoothing threshold
  big_df <- data.frame(
    position = seq_len(6000),
    mean     = rnorm(6000, 18, 2),
    median   = rnorm(6000, 17, 2),
    q25      = rnorm(6000, 14, 2),
    q75      = rnorm(6000, 21, 2)
  )
  p <- make_per_position_plot(big_df, file_name = "test.fastq")
  expect_s3_class(p, "ggplot")
})

test_that("make_per_position_plot uses file_name as plot title", {
  p <- make_per_position_plot(make_fake_per_pos(), file_name = "my_sample.fastq")
  expect_equal(p$labels$title, "my_sample.fastq")
})
