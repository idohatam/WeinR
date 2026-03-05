library(testthat)
library(WeinR)

# Minimal fake object for tests
fake_mfa <- function() {
  structure(list(dummy = TRUE), class = "not_LongReadQC")
}

test_that("Reporter errors when mfa is NULL", {
  
  expect_error(
    Reporter(
      path = tempfile(fileext = ".Rmd"),
      render_html = FALSE,
      mfa = NULL
    )
  )
  
})

test_that("Reporter writes an Rmd file and returns path", {
  
  tmp <- tempdir()
  out <- file.path(tmp, "report")
  
  res <- Reporter(
    path = out,
    render_html = FALSE,
    overwrite = TRUE,
    mfa = fake_mfa()
  )
  
  expect_true(file.exists(paste0(out, ".Rmd")))
  expect_type(res, "list")
  expect_true("rmd" %in% names(res))
})

test_that("Reporter appends .Rmd when extension missing", {
  
  tmp <- tempdir()
  out <- file.path(tmp, "my_report")
  
  Reporter(
    path = out,
    render_html = FALSE,
    overwrite = TRUE,
    mfa = fake_mfa()
  )
  
  expect_true(file.exists(paste0(out, ".Rmd")))
})

test_that("Reporter prevents overwrite when overwrite = FALSE", {
  
  tmp <- tempdir()
  out <- file.path(tmp, "report.Rmd")
  
  writeLines("existing file", out)
  
  expect_error(
    Reporter(
      path = out,
      render_html = FALSE,
      overwrite = FALSE,
      mfa = fake_mfa()
    ),
    "File already exists"
  )
  
})

test_that("Reporter writes expected YAML fields", {
  
  tmp <- tempdir()
  out <- file.path(tmp, "report.Rmd")
  
  Reporter(
    path = out,
    title = "QC Title",
    author = "Owen",
    date = "March 4, 2026",
    theme = "cosmo",
    highlight = "tango",
    code_folding = "hide",
    render_html = FALSE,
    overwrite = TRUE,
    mfa = fake_mfa()
  )
  
  txt <- paste(readLines(out, warn = FALSE), collapse = "\n")
  
  expect_match(txt, 'title: "QC Title"', fixed = TRUE)
  expect_match(txt, 'author: "Owen"', fixed = TRUE)
  expect_match(txt, 'date: "March 4, 2026"', fixed = TRUE)
  expect_match(txt, "theme: cosmo", fixed = TRUE)
  expect_match(txt, "highlight: tango", fixed = TRUE)
  expect_match(txt, "code_folding: hide", fixed = TRUE)
  
})

test_that("Reporter validates code_folding argument", {
  
  expect_error(
    Reporter(
      path = tempfile(fileext = ".Rmd"),
      code_folding = "bad_option",
      render_html = FALSE,
      mfa = fake_mfa()
    )
  )
  
})

test_that("Reporter includes metadata section only when metadata = TRUE", {
  
  tmp <- tempdir()
  
  out1 <- file.path(tmp, "no_metadata.Rmd")
  
  Reporter(
    path = out1,
    render_html = FALSE,
    overwrite = TRUE,
    metadata = FALSE,
    mfa = fake_mfa()
  )
  
  txt1 <- paste(readLines(out1, warn = FALSE), collapse = "\n")
  
  expect_false(grepl("## Metadata", txt1, fixed = TRUE))
  
  out2 <- file.path(tmp, "with_metadata.Rmd")
  
  Reporter(
    path = out2,
    render_html = FALSE,
    overwrite = TRUE,
    metadata = TRUE,
    mfa = fake_mfa()
  )
  
  txt2 <- paste(readLines(out2, warn = FALSE), collapse = "\n")
  
  expect_true(grepl("## Metadata", txt2, fixed = TRUE))
  
})

test_that("Reporter errors if CSS asset is missing", {
  
  local_mocked_bindings(
    system.file = function(...) ""
  )
  
  expect_error(
    Reporter(
      path = tempfile(fileext = ".Rmd"),
      render_html = FALSE,
      overwrite = TRUE,
      mfa = fake_mfa()
    ),
    "Could not find 'styles/colorblind.css'"
  )
  
})

test_that("Reporter errors if legend PNG is missing", {
  
  calls <- 0L
  
  local_mocked_bindings(
    system.file = function(...) {
      
      calls <<- calls + 1L
      
      # First call = CSS (valid)
      if (calls == 1) return("/fake/css.css")
      
      # Second call = missing PNG
      ""
      
    }
  )
  
  expect_error(
    Reporter(
      path = tempfile(fileext = ".Rmd"),
      render_html = FALSE,
      overwrite = TRUE,
      mfa = fake_mfa()
    ),
    "quality_legend"
  )
  
})