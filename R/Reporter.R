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
  
  # Locate package CSS (installed from inst/styles/colorblind.css)
  pkg_css <- system.file("styles", "colorblind.css", package = "WeinR")
  if (!nzchar(pkg_css)) {
    stop("Could not find 'styles/colorblind.css' in the installed WeinR package.", call. = FALSE)
  }
  
  # DEV MODE css (absolute path is fine)
  #pkg_css <- normalizePath(
   # file.path("inst", "styles", "colorblind.css"),
    #winslash = "/",
    #mustWork = TRUE
  #)
  
  
  # Legend PNGs (absolute paths)
  quality_legend_png <- system.file("images", "quality_legend.png", package = "WeinR")
  if (!nzchar(pkg_css)) {
    stop("Could not find 'images/quality_legend.png' in the installed WeinR package.", call. = FALSE)
  }
  #quality_legend_png <- normalizePath(
  # file.path("inst", "images", "quality_legend.png"),
  #winslash = "/",
  #mustWork = TRUE
  #)
  
  per_pos_q_legend_png <- system.file("images", "per_pos_q_legend.png", package = "WeinR")
  if (!nzchar(pkg_css)) {
    stop("Could not find 'images/per_pos_q_legend.png' in the installed WeinR package.", call. = FALSE)
  }
  #per_pos_q_legend_png <- normalizePath(
  # file.path("inst", "images", "per_pos_q_legend.png"),
  #  winslash = "/",
  #  mustWork = TRUE
  #)
  
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
    "stopifnot(requireNamespace('ggplot2', quietly = TRUE))",
    "stopifnot(requireNamespace('DT', quietly = TRUE))",
    "stopifnot(requireNamespace('gridExtra', quietly = TRUE))",
    "stopifnot(requireNamespace('grid', quietly = TRUE))",
    "",
    "# Render a grid of ggplots OR image-file paths to a single PNG, then include it",
    "render_grid_png <- function(grobs, max_cols = 3, height_in = 3.5, dpi = 300) {",
    "  if (is.null(grobs)) return(invisible(NULL))",
    "",
    "  # Flatten nested lists",
    "  if (is.list(grobs) && !inherits(grobs, 'gg')) {",
    "    grobs_flat <- unlist(grobs, recursive = TRUE, use.names = FALSE)",
    "  } else {",
    "    grobs_flat <- grobs",
    "  }",
    "",
    "  # If these are file paths -> build rasterGrobs, then arrange",
    "  if (is.character(grobs_flat)) {",
    "    paths <- grobs_flat[!is.na(grobs_flat) & nzchar(grobs_flat) & file.exists(grobs_flat)]",
    "    if (!length(paths)) return(invisible(NULL))",
    "",
    "    stopifnot(requireNamespace('png', quietly = TRUE))",
    "",
    "    rasters <- lapply(paths, function(p) {",
    "      img <- png::readPNG(p)",
    "      grid::rasterGrob(img, interpolate = TRUE)",
    "    })",
    "",
    "    n <- length(rasters)",
    "    ncol <- min(max_cols, n)",
    "    nrow <- ceiling(n / ncol)",
    "    width_in <- if (ncol == 1) 8 else if (ncol == 2) 11 else 14",
    "",
    "    tmp <- tempfile(fileext = '.png')",
    "    grDevices::png(filename = tmp, width = width_in * dpi, height = height_in * nrow * dpi, res = dpi)",
    "    gridExtra::grid.arrange(grobs = rasters, ncol = ncol, nrow = nrow)",
    "    grDevices::dev.off()",
    "    return(knitr::include_graphics(tmp))",
    "  }",
    "",
    "  # Otherwise assume plot grobs/ggplots",
    "  n <- length(grobs)",
    "  if (!n) return(invisible(NULL))",
    "  ncol <- min(max_cols, n)",
    "  nrow <- ceiling(n / ncol)",
    "  width_in <- if (ncol == 1) 8 else if (ncol == 2) 11 else 14",
    "",
    "  tmp <- tempfile(fileext = '.png')",
    "  grDevices::png(filename = tmp, width = width_in * dpi, height = height_in * nrow * dpi, res = dpi)",
    "  gridExtra::grid.arrange(grobs = grobs, ncol = ncol, nrow = nrow)",
    "  grDevices::dev.off()",
    "  knitr::include_graphics(tmp)",
    "}",
    "```",
    "",
    "<div class='section-card'>",
    "## **Summary Statistics by File**",
    "```{r}",
    "summ <- qc@summary_metrics",
    "if (!is.null(summ) && NROW(summ) > 0 && 'file' %in% names(summ)) {",
    "  summ <- summ[!is.na(summ$file) & nzchar(trimws(summ$file)), , drop = FALSE]",
    "}",
    "if (!is.null(summ) && NROW(summ) > 0) {",
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
    "} else { cat('No summary metrics available.') }",
    "```",
    "</div>",
    "",
    "<hr class='section-sep'/>",
    ""
  )
  
  if (isTRUE(metadata)) {
    body <- c(
      body,
      "<div class='section-card'>",
      "## Metadata",
      "```{r}",
      "md <- qc@metadata",
      "",
      "get_step_tbl <- function(md, step_name, pretty_step) {",
      "  files <- names(md)",
      "  if (!length(files)) return(NULL)",
      "",
      "  reads_before <- vapply(files, function(f) {",
      "    s <- md[[f]][[step_name]]",
      "    if (is.null(s) || is.null(s[['reads_before']])) return(NA_integer_)",
      "    as.integer(s[['reads_before']])",
      "  }, integer(1))",
      "",
      "  reads_after <- vapply(files, function(f) {",
      "    s <- md[[f]][[step_name]]",
      "    if (is.null(s) || is.null(s[['reads_after']])) return(NA_integer_)",
      "    as.integer(s[['reads_after']])",
      "  }, integer(1))",
      "",
      "  data.frame(",
      "    Filename = files,",
      "    Step = pretty_step,",
      "    `Reads Before` = reads_before,",
      "    `Reads After` = reads_after,",
      "    check.names = FALSE,",
      "    stringsAsFactors = FALSE",
      "  )",
      "}",
      "",
      "if (length(md)) {",
      "  tbl_filter  <- get_step_tbl(md, 'filter_summary',  'Filter')",
      "  tbl_adapter <- get_step_tbl(md, 'adapter_summary', 'Adapter')",
      "  tbl_trim    <- get_step_tbl(md, 'trim_summary',    'Trim')",
      "  meta_tbl <- do.call(rbind, Filter(Negate(is.null), list(tbl_filter, tbl_adapter, tbl_trim)))",
      "  meta_tbl <- meta_tbl[, c('Filename','Step','Reads Before','Reads After'), drop = FALSE]",
      "  DT::datatable(",
      "    meta_tbl,",
      "    options = list(",
      "      scrollX = TRUE,",
      "      pageLength = 25,",
      "      autoWidth = TRUE,",
      "      order = list(list(0, 'asc'))",
      "    ),",
      "    rownames = FALSE,",
      "    caption = 'Processing metadata',",
      "    class = 'stripe hover'",
      "  ) |>",
      "  DT::formatRound(columns = c('Reads Before','Reads After'), digits = 0, mark = ',')",
      "} else { cat('No metadata available.') }",
      "```",
      "</div>",
      "",
      "<hr class='section-sep'/>",
      ""
    )
  }
  
  body <- c(
    body,
    "<div class='section-card'>",
    "## **Plots**",
    "",
    "### **Read Length Distribution**",
    "```{r, fig.align='left', out.width='70%'}",
    "lens <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['length_hist']] else NULL))",
    "if (length(lens)) render_grid_png(lens, max_cols = 3, height_in = 3.5) else cat('No length histograms available')",
    "```",
    "",
    "### **Average Q-score Distribution**",
    "```{r, fig.align='left', out.width='20%'}",
    sprintf("knitr::include_graphics('%s')", quality_legend_png),
    "```",
    "```{r, fig.align='left', out.width='70%'}",
    "qhist <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['quality_hist']] else NULL))",
    "if (length(qhist)) render_grid_png(qhist, max_cols = 3, height_in = 3.5) else cat('No Q-score histograms available')",
    "```",
    "",
    "### **GC-content Distribution**",
    "```{r, fig.align='left', out.width='70%'}",
    "gch <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['gc_hist']] else NULL))",
    "if (length(gch)) render_grid_png(gch, max_cols = 3, height_in = 3.5) else cat('No GC-content histograms available')",
    "```",
    "",
    "### **Read Quality vs Read Length**",
    "```{r, fig.align='left', out.width='70%'}",
    "qvl <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['q_vs_length']] else NULL))",
    "if (length(qvl)) render_grid_png(qvl, max_cols = 3, height_in = 3.5) else cat('No quality-vs-length plots available')",
    "```",
    "",
    "### **Per-position quality (binned at 2000 bp)**",
    "```{r, fig.align='left', out.width='20%'}",
    sprintf("knitr::include_graphics('%s')", per_pos_q_legend_png),
    "```",
    "```{r, fig.align='left', out.width='70%'}",
    "ppq <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['per_pos_q']] else NULL))",
    "if (length(ppq)) render_grid_png(ppq, max_cols = 3, height_in = 3.5) else cat('No per-position quality plots available')",
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
      stop("Package 'rmarkdown' is required to render HTML. Install it via install.packages('rmarkdown').",
           call. = FALSE)
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
