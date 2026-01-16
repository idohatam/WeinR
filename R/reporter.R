#' Create an HTML QC report from a LongReadQC object
#' Writes an R Markdown (`.Rmd`) file to `path` and (optionally) renders it to
#' an HTML report using `rmarkdown::render()`. The report includes a DT summary
#' table and plot sections assembled from `mfa` (a `LongReadQC` object).
#'
#' @param path Character. Output path for the report **without or with** the `.Rmd`
#'   extension. If no extension is provided, `.Rmd` is appended. Directories are
#'   created if needed.
#' @param title Character. Report title shown in the rendered HTML.
#' @param author Character. Report author shown in the rendered HTML.
#' @param date Character. Report date shown in the rendered HTML.
#' @param toc Logical. Whether to include a table of contents.
#' @param code_folding Character. One of `"none"`, `"show"`, `"hide"`. Controls
#'   code folding in the HTML output.
#' @param theme Character. R Markdown Bootstrap theme name (e.g. `"cosmo"`).
#' @param highlight Character. Syntax highlighting theme (e.g. `"tango"`).
#' @param overwrite Logical. If `FALSE` and the output `.Rmd` exists, an error is
#'   thrown. If `TRUE`, existing files are replaced.
#' @param render_html Logical. If `TRUE`, renders the `.Rmd` to HTML via
#'   `rmarkdown::render()`.
#' @param open_browser Logical. If `TRUE` and `render_html=TRUE`, opens the
#'   rendered HTML in a browser.
#' @param mfa A `LongReadQC` object (required). Used by the report as
#'   `params$mfa`.
#'
#' @return If `render_html=TRUE`, an (invisible) list with elements `rmd` and
#'   `html` giving normalized paths to the generated files. If `render_html=FALSE`,
#'   an (invisible) list with element `rmd`.
#'
#' @details
#' This function assumes a stylesheet exists at `<project root>/styles/colorblind.css`,
#' where `<project root>` is the current working directory (`getwd()`).
#'
#' The generated report expects `mfa` to be a `LongReadQC` S4 object with fields
#' used in the report body (e.g., `@summary_metrics`, `@plots`, `@files`).

CreateReport <- function(
    path          = "report.Rmd",
    title         = "My Report",
    author        = Sys.info()[["user"]],
    date          = format(Sys.Date(), "%B %d, %Y"),
    toc           = TRUE,
    code_folding  = c("none", "show", "hide"),
    theme         = "cosmo",
    highlight     = "tango",
    overwrite     = FALSE,
    render_html   = TRUE,
    open_browser  = interactive(),
    mfa           = NULL
) {
  stopifnot(!is.null(mfa))
  code_folding <- match.arg(code_folding)
  
  # normalize and ensure .Rmd suffix
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!grepl("\\.Rmd$", path, ignore.case = TRUE)) {
    path <- paste0(path, ".Rmd")
  }
  
  # ensure output directory exists
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  # overwrite guard
  if (file.exists(path) && !overwrite) {
    stop(
      "File already exists: ", normalizePath(path),
      "\nSet overwrite = TRUE to replace it."
    )
  }
  
  css_path <- normalizePath(
    file.path(getwd(), "styles", "colorblind.css"),
    winslash = "/",
    mustWork = FALSE
  )
  
  yaml <- c(
    "---",
    sprintf('title: "%s"', title),
    sprintf('author: "%s"', author),
    sprintf('date: "%s"', date),
    "output:",
    "  html_document:",
    sprintf("    toc: %s", tolower(as.character(toc))),
    "    toc_depth: 3",
    sprintf("    theme: %s", theme),
    sprintf("    highlight: %s", highlight),
    sprintf("    code_folding: %s", code_folding),
    "    df_print: paged",
    "    self_contained: true",
    sprintf('    css: "%s"', css_path),
    "params:",
    "  mfa: !r NULL",
    "---",
    ""
  )
  
  body <- c(
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
    "knitr::opts_knit$set(root.dir = dirname(knitr::current_input()))",
    "qc <- params$mfa",
    "stopifnot(methods::is(qc, 'LongReadQC'))",
    "suppressPackageStartupMessages({",
    "  library(ggplot2)",
    "  library(DT)",
    "  library(gridExtra)",
    "})",
    "",
    "# Render a gridExtra arrangement to a PNG at a chosen width, then include it",
    "render_grid_png <- function(grobs, max_cols = 3, height_in = 4, dpi = 120) {",
    "  n <- length(grobs)",
    "  if (!n) return(invisible(NULL))",
    "  ncol <- min(max_cols, n)",
    "  nrow <- ceiling(n / ncol)",
    "",
    "  width_in <- if (ncol == 1) 10 else if (ncol == 2) 12 else 30",
    "",
    "  tmp <- tempfile(fileext = '.png')",
    "  grDevices::png(filename = tmp, width = width_in * dpi, height = height_in * dpi, res = dpi)",
    "  gridExtra::grid.arrange(grobs = grobs, ncol = ncol, nrow = nrow)",
    "  grDevices::dev.off()",
    "",
    "  knitr::include_graphics(tmp)",
    "}",
    "```",
    "",
    "<div class='section-card'>",
    "## **Summary Statistics by File**",
    "```{r}",
    "summ <- qc@summary_metrics",
    "",
    "# Drop rows with missing or empty file names",
    "if (!is.null(summ) && NROW(summ) > 0 && 'file' %in% names(summ)) {",
    "  summ <- summ[!is.na(summ$file) & nzchar(trimws(summ$file)), , drop = FALSE]",
    "}",
    "",
    "if (!is.null(summ) && NROW(summ) > 0) {",
    "  # Ensure file column is first",
    "  if ('file' %in% names(summ)) {",
    "    cols <- c('file', setdiff(names(summ), 'file'))",
    "    summ <- summ[, cols, drop = FALSE]",
    "  }",
    "  DT::datatable(",
    "    summ,",
    "    options = list(scrollX = TRUE, autoWidth = TRUE),",
    "    rownames = FALSE,",
    "    caption = 'Per-file summary metrics',",
    "    class = 'stripe hover'",
    "  )",
    "} else {",
    "  cat('No summary metrics available.')",
    "}",
    "```",
    "</div>",
    "",
    "<hr class='section-sep'/>",
    "",
    "<div class='section-card'>",
    "## **Plots**",
    "",
    "### **Read Length Distribution**",
    "```{r, fig.align='center', out.width='100%'}",
    "lens <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['length_hist']] else NULL))",
    "if (length(lens)) {",
    "  render_grid_png(lens, max_cols = 3, height_in = 4)",
    "} else { cat('No length histograms available') }",
    "```",
    "",
    "### **Average Q-score Distribution**",
    "```{r, fig.align='center', out.width='100%'}",
    "qhist <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['quality_hist']] else NULL))",
    "if (length(qhist)) {",
    "  render_grid_png(qhist, max_cols = 3, height_in = 4)",
    "} else { cat('No Q-score histograms available') }",
    "```",
    "",
    "### **GC-content Distribution**",
    "```{r, fig.align='center', out.width='100%'}",
    "gch <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['gc_hist']] else NULL))",
    "if (length(gch)) {",
    "  render_grid_png(gch, max_cols = 3, height_in = 4)",
    "} else { cat('No GC-content histograms available') }",
    "```",
    "",
    "### **Read Quality vs Read Length**",
    "```{r, fig.align='center', out.width='100%'}",
    "qvl <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['q_vs_length']] else NULL))",
    "if (length(qvl)) {",
    "  render_grid_png(qvl, max_cols = 3, height_in = 4)",
    "} else { cat('No quality-vs-length plots available') }",
    "```",
    "",
    "### **Per-position quality  Binned at 2000 bp**",
    "```{r, fig.align='center', out.width='100%'}",
    "ppq <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['per_pos_q']] else NULL))",
    "if (length(ppq)) {",
    "  render_grid_png(ppq, max_cols = 3, height_in = 4)",
    "} else { cat('No per-position quality plots available') }",
    "```",
    "</div>",
    "",
    "<hr class='section-sep'/>",
    "",
    "<div class='section-card'>",
    "## **Session info**",
    "```{r}",
    "sessionInfo()",
    "```",
    "</div>",
    "",
    "<div class='page-footer'>",
    "  Generated With WeinR",
    "</div>"
  )
  
  writeLines(c(yaml, body), con = path)
  
  if (isTRUE(render_html)) {
    if (!requireNamespace('rmarkdown', quietly = TRUE)) {
      stop("Package 'rmarkdown' is required to render HTML. Install it via install.packages('rmarkdown').")
    }
    
    html_out <- rmarkdown::render(
      path,
      output_format = "html_document",
      params        = list(mfa = mfa),
      knit_root_dir = dirname(path),
      envir         = new.env(parent = globalenv()),
      quiet         = TRUE
    )
    
    if (open_browser) utils::browseURL(html_out)
    return(invisible(list(rmd = normalizePath(path), html = normalizePath(html_out))))
  }
  
  invisible(list(rmd = normalizePath(path)))
}
