#' Create an HTML QC report from a LongReadQC object
#'
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
#' @param metadata Logical. If `TRUE`, include a metadata table summarizing
#'
#' @return If `render_html=TRUE`, an (invisible) list with elements `rmd` and
#'   `html` giving normalized paths to the generated files. If `render_html=FALSE`,
#'   an (invisible) list with element `rmd`.
#'
#' @details
#' The stylesheet is shipped with the package at `inst/styles/colorblind.css`.
#' For portability, it is copied next to the generated `.Rmd` and referenced
#' using a relative path in the YAML header.
Reporter <- function(
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
    mfa           = NULL,
    metadata      = FALSE
) {
  stopifnot(!is.null(mfa))
  code_folding <- match.arg(code_folding)
  
  # Normalize and ensure .Rmd suffix
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!grepl("\\.Rmd$", path, ignore.case = TRUE)) {
    path <- paste0(path, ".Rmd")
  }
  
  # Ensure output directory exists
  out_dir <- dirname(path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Overwrite guard
  if (file.exists(path) && !overwrite) {
    stop(
      "File already exists: ", normalizePath(path),
      "\nSet overwrite = TRUE to replace it.",
      call. = FALSE
    )
  }
  ###########
  ###########
  ###########
  ###########
  # Locate package CSS (installed from inst/styles/colorblind.css)
  #pkg_css <- system.file("styles", "colorblind.css", package = "WeinR")
  #if (!nzchar(pkg_css)) {
  #  stop("Could not find 'styles/colorblind.css' in the installed WeinR package.", call. = FALSE)
  #}
  ###########
  ###########
  ###########
  ###########
  # ---- DEV MODE: load CSS directly from inst/ ----
  # NOTE: This is for development only. Replace with system.file()
  # before release / Bioconductor submission.
  pkg_css <- normalizePath(file.path("inst", "styles", "colorblind.css"),
                           winslash = "/", mustWork = TRUE)
  
  yaml <- c(
    "---",
    sprintf('title: "%s"', title),
    sprintf('author: "%s"', author),
    sprintf('date: "%s"', date),
    "output:",
    "  html_document:",
    sprintf("    theme: %s", theme),
    sprintf("    highlight: %s", highlight),
    sprintf("    code_folding: %s", code_folding),
    "    df_print: paged",
    "    self_contained: true",
    #    '    css: "colorblind.css"',
    sprintf('    css: "%s"', pkg_css),
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
    "",
    "# Require suggested packages only when rendering the report",
    "stopifnot(requireNamespace('ggplot2', quietly = TRUE))",
    "stopifnot(requireNamespace('DT', quietly = TRUE))",
    "stopifnot(requireNamespace('gridExtra', quietly = TRUE))",
    "",
    "# Render a gridExtra arrangement to a PNG at a chosen width, then include it",
    "render_grid_png <- function(grobs, max_cols = 3, height_in = 4, dpi = 120) {",
    "  n <- length(grobs)",
    "  if (!n) return(invisible(NULL))",
    "  # ---- NEW: if grobs are file paths, include directly ----",
    "  if (is.list(grobs)) {",
    "    # common case: list of character paths (or nested lists of character paths)",
    "    flat <- unlist(grobs, use.names = FALSE)",
    "    if (is.character(flat)) {",
    "      flat <- flat[!is.na(flat) & nzchar(flat) & file.exists(flat)]",
    "      if (!length(flat)) return(invisible(NULL))",
    "      return(knitr::include_graphics(flat))",
    "    }",
    "  }",
    "  if (is.character(grobs)) {",
    "    grobs <- grobs[!is.na(grobs) & nzchar(grobs) & file.exists(grobs)]",
    "    if (!length(grobs)) return(invisible(NULL))",
    "    return(knitr::include_graphics(grobs))",
    "  }",
    "",
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
    ""
  )
  if (metadata == TRUE) {
    body <- c(
      body,
      "<div class='section-card'>",
      "## Metadata",
      "```{r}",
      "md <- qc@metadata",
      "",
      "if (length(md)) {",
      "  files <- names(md)",
      "",
      "  reads_before <- vapply(files, function(f) {",
      "    fs <- md[[f]][['filter_summary']]",
      "    if (is.null(fs) || is.null(fs[['reads_before']])) return(NA_integer_)",
      "    as.integer(fs[['reads_before']])",
      "  }, integer(1))",
      "",
      "  reads_after <- vapply(files, function(f) {",
      "    fs <- md[[f]][['filter_summary']]",
      "    if (is.null(fs) || is.null(fs[['reads_after']])) return(NA_integer_)",
      "    as.integer(fs[['reads_after']])",
      "  }, integer(1))",
      "",
      "  filtered_fname <- vapply(files, function(f) {",
      "    fs <- md[[f]][['filter_summary']]",
      "    op <- if (!is.null(fs)) fs[['output_paths']] else NULL",
      "    if (is.null(op) || length(op) == 0 || is.na(op[1])) return(NA_character_)",
      "    basename(op[1])  # substring after last '/'",
      "  }, character(1))",
      "",
      "  meta_tbl <- data.frame(",
      "    Filename = files,",
      "    `Pre-filter Reads` = reads_before,",
      "    `Post-filter Reads` = reads_after,",
      "    `Filtered Filename` = filtered_fname,",
      "    check.names = FALSE,",
      "    stringsAsFactors = FALSE",
      "  )",
      "",
      "  DT::datatable(",
      "    meta_tbl,",
      "    options = list(scrollX = TRUE, pageLength = 15, autoWidth = TRUE),",
      "    rownames = FALSE,",
      "    caption = 'Filtering metadata',",
      "    class = 'stripe hover'",
      "  ) |>",
      "  DT::formatRound(",
      "    columns = c('Pre-filter Reads', 'Post-filter Reads'),",
      "    digits = 0,",
      "    mark = ','",
      "  )",
      "} else {",
      "  cat('No metadata available.')",
      "}",
      "```",
      "</div>"
    )
  }
  body <- c(
    body,
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
    "### **Per-position quality (binned at 2000 bp)**",
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
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
      stop(
        "Package 'rmarkdown' is required to render HTML. Install it via install.packages('rmarkdown').",
        call. = FALSE
      )
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
    return(invisible(list(
      rmd  = normalizePath(path, winslash = "/", mustWork = FALSE),
      html = normalizePath(html_out, winslash = "/", mustWork = FALSE)
    )))
  }
  
  invisible(list(rmd = normalizePath(path, winslash = "/", mustWork = FALSE)))
}
